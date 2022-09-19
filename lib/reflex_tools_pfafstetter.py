"""
Library Features:

Name:          reflex_tools_pfafstetter
Author(s):     Mauro Arcorace (mauro.arcorace@cimafoundation.org)
               Alessandro Masoero (alessandro.masoero@cimafoundation.org)
               Valerio Basso
               Giulia Bruno (giulia.bruno@cimafoundation.org)
               Alessia Matanò
               Andrea Libertino (andrea.libertino@cimafoundation.org)
Date:          '20220916'
Version:       '2.0.1'
"""
########################################################################################################################
# Import libraries
import fiona
import fiona.crs
import geopandas as gpd
import json
import numpy as np
import os
from osgeo import osr
import pandas as pd
import subprocess
import sys
from shapely.geometry import shape
########################################################################################################################

# Function for assign pfafstetter code
def assign_pfafstetter_code (sPathInput, sDomain, input_epsg, in_filename, iNum_col, sPathOutput):

    iNum_col = int(iNum_col)
    # path
    # sPathInput = sPath
    # sPathInput = os.path.join(sPath, in_dirname)
    # sPathInput = os.path.join(sPath, 'input')
    # sPathOutput = os.path.join(sPath, out_dirname)

    # read gpkg
    # sInput = os.path.join(sPathInput, sDomain + '_str_order_90m.gpkg')
    sInput = os.path.join(sPathInput, in_filename)

    # sInputMB = os.path.join(sPathInput, sDomain + '_mbas_90m.gpkg')

    pd.set_option('mode.chained_assignment', None)

    collection = list(fiona.open(sInput, 'r'))
    df1 = pd.DataFrame(collection)
    # df1['isline'] = np.nan
    # for i in np.arange(0, len(df1), 1):
    #    if df1['geometry'][i]['type'] == 'Point':
    #        df1['isline'][i] = 0
    #    else:
    #        df1['isline'][i] = 1
    # df1 = df1[df1['isline'] == 1]  # read just line string geometry

    df1['isvalid'] = df1['geometry'].apply(lambda x: isvalid(x))
    idx = df1.index[df1['isvalid'] == 0]
    for i in np.arange(0, len(idx), 1):
        a = df1['geometry'][idx[i]]['coordinates']
        df1['geometry'][idx[i]]['coordinates'] = [a[0], a[0]]
    collection = json.loads(df1.to_json(orient='records'))
    gdf_orig = gpd.GeoDataFrame.from_features(collection)  # Convert to geodataframe
    # gdf_orig = gdf_orig.drop(columns='cat')
    gdf_orig = gdf_orig.rename(columns={"cat": "cat_bas"})
    gdf_orig.columns = gdf_orig.columns.str.replace('a_', '')  # NEW line to remove 'a_' prefix in columns header
    # gdf_orig = [col.replace('a_', '') for col in gdf_orig]##########-----------------------------NEW line
    # gdf_STR = gpd.GeoDataFrame.from_features(collection)  # Convert to geodataframe
    # gdf_STR.crs = {'init' :'epsg:4326'}
    # gdf_MB = gpd.read_file(sInputMB)
    # gdf_MB.crs = {'init' :'epsg:4326'}
    # gdf_orig = gpd.sjoin(gdf_STR, gdf_MB, op='intersects')

    # Pfafstetter labelling
    gdf_orig = gdf_orig.assign(stream_d=gdf_orig['stream'])
    gdf_orig = gdf_orig.dissolve(by='stream_d')
    gdf_orig.loc[~gdf_orig['next_strea'].isin(gdf_orig['stream']), 'next_strea'] = -1  # set as outlet streams whose next streams are not in the gpkg (points and not lines)
    gdf_orig['pfA'] = np.nan  # add blank column to geodataframe
    gdf_orig['FiscalCode'] = np.nan
    gdf_orig['MBID_new'] = np.nan
    gdf_bkp = gdf_orig.copy(deep=True)  # create a backup dataframe, used at the end of the script
    cc = gdf_orig['MBID'].unique()
    ccs = cc.size
    count = 0

    for m in np.arange(0, ccs, 1):  # for loop over each macro basin
        count = count + 1
        gdf_orig_mb = gdf_orig.loc[gdf_orig['MBID'] == cc[m]]
        gdf_orig_mb = gdf_orig_mb.drop_duplicates(subset='stream')  # delete duplicated and null streams
        gdf_orig_mb = gdf_orig_mb.dropna(subset=['stream'])
        # print('working on macrobasin ' + str(cc[m]) + ' (' + str(count) + '/' + str(ccs) + ')')

        outlet = gdf_orig_mb.loc[gdf_orig_mb['next_strea'] == -1, 'stream']  # find outlet stream
        gdf_outlet = gdf_orig_mb.loc[gdf_orig_mb['next_strea'] == -1]
        gdf_no_outlet = gdf_orig_mb.loc[gdf_orig_mb['next_strea'] != -1]

        for y in np.arange(0, gdf_outlet['stream'].unique().size, 1):  # loop over outlets
            # print('working on outlet ' + str(gdf_outlet['stream'].unique()[y]))
            gdf_to_add = gdf_outlet.loc[gdf_outlet['stream'] == gdf_outlet['stream'].unique()[y]]
            gdf_final = gdf_to_add.append(gdf_no_outlet)

            tovisit = []
            iHack_max = gdf_final['hack'].max()
            # print('Hack max = ' + str(iHack_max))

            for j in np.arange(1, iHack_max + 1, 1):  # for loop over needed digits to have unique code
                # print('ITERATION ' + str(j))
                if j == 1:
                    gdf = gdf_final
                    outlet_final = gdf_final.loc[gdf_final['next_strea'] == -1, 'stream']
                    outlet_final = outlet_final.reset_index(drop=True)
                    dbjun = get_junction_streams(gdf)  # extract junction streams
                    iHack_outlet = gdf.loc[gdf['next_strea'] == -1, 'hack'].squeeze()  # extract main tributaries
                    main_str = select_mainstreams(dbjun, iHack_outlet + 1)

                    tovisit = outlet_final  # initialise
                    gdf.loc[gdf['stream'] == outlet_final[0], 'pfA'] = 1* 10 ** (iHack_max - j)
                    i = 1
                    while 1 <= tovisit.size:  # while loop
                        visiting = tovisit.values.astype(int)[0]
                        n = explore_tree(gdf, visiting, iNum_col)  # streams to explore
                        # print(n.values)
                        if sum(n) != 0:
                            if n[n.shape[0] - 1] == 0:
                                n = n[0:n.shape[0] - 1]
                            else:
                                # print('previous streams exist')
                                pass
                            a = pd.Series(list(
                                set(n) & set(main_str)))  # intersection between streams to explore and main streams
                            if a.size == 0:  # normal junction (same Pfafstetter index)
                                gdf.loc[gdf['stream'].isin(n), 'pfA'] = (gdf.loc[
                                    gdf['stream'] == visiting, 'pfA'].values.astype(int)[0])
                                tovisit = tovisit.append(n, ignore_index=True)
                                tovisit_temp = tovisit.to_list()
                                tovisit_temp.pop(0)
                                tovisit = pd.Series(tovisit_temp)
                            else:  # main junction (main tributary gets +1 [even number], main channel gets +2 [odd number])
                                gdf.loc[gdf['stream'].isin(n), 'pfA'] = (gdf.loc[
                                    gdf['stream'] == visiting, 'pfA'].values.astype(int)[0]) + 2 * 10 ** (iHack_max - j)
                                gdf.loc[gdf['stream'].isin(a), 'pfA'] = (gdf.loc[
                                    gdf['stream'] == visiting, 'pfA'].values.astype(int)[0]) + 1 * 10 ** (iHack_max - j)
                                tovisit = tovisit.append(n, ignore_index=True)
                                tovisit_temp = tovisit.to_list()
                                tovisit_temp.pop(0)
                                tovisit = pd.Series(tovisit_temp)
                        else:  # source stream, go back to junction
                            tovisit_temp = tovisit.to_list()
                            tovisit_temp.pop(0)
                            tovisit = pd.Series(tovisit_temp)
                        i = i + 1  # number of iterations
                    for k in np.arange(0, gdf['pfA'].unique().size,
                                       1):  # loop over subbasins (to explore in the next step)
                        gdf_tmp = gdf.loc[gdf['pfA'] == gdf['pfA'].unique()[k]]
                        gdf_tmp.loc[~gdf_tmp['next_strea'].isin(gdf_tmp['stream']), 'next_strea'] = -1  # outlet
                        for i in np.arange(1, gdf.shape[1] - iNum_col):
                            gdf_tmp.loc[~gdf_tmp['prev_str0' + str(i)].isin(gdf_tmp['stream']), 'prev_str0' + str(
                                i)] = 0  # source stream
                        gdf.loc[gdf['pfA'] == gdf['pfA'].unique()[k]] = gdf_tmp
                    gdf_final = gdf

                else:
                    aa = gdf_final['pfA'].unique()
                    for k in np.arange(0, aa.size, 1):
                        gdf = gdf_final.loc[gdf_final['pfA'] == aa[k]]
                        if len(gdf.index) <= 1:
                            pass
                        else:
                            outlet = gdf.loc[gdf['next_strea'] == -1, 'stream']  # find outlet stream
                            dbjun = get_junction_streams(gdf)  # extract junction streams
                            iHack_outlet = gdf.loc[gdf['next_strea'] == -1, 'hack'].values.astype(int)[
                                0]  # extract main tributaries
                            main_str = select_mainstreams(dbjun, iHack_outlet + 1)
                            tovisit = outlet  # initialise
                            gdf.loc[gdf['stream'] == outlet.values.astype(int)[0], 'pfA'] = gdf.loc[gdf['stream'] ==
                                                                                                    outlet.values.astype(
                                                                                                        int)[
                                                                                                        0], 'pfA'] + 1 * 10 ** (
                                                                                                        iHack_max - j)
                            i = 1
                            while 1 <= tovisit.size:
                                visiting = tovisit.values.astype(int)[0]
                                n = explore_tree(gdf, visiting, iNum_col)  # streams to explore
                                # print(n.values)
                                if sum(n) != 0:
                                    if n[n.shape[0] - 1] == 0:
                                        n = n[0:n.shape[0] - 1]
                                    else:
                                        pass
                                    # print('previous streams exist')
                                    a = pd.Series(list(set(n) & set(
                                        main_str)))  # intersection between streams to explore and main streams
                                    if a.size == 0:  # normal junction (same Pfafstetter index)
                                        gdf.loc[gdf['stream'].isin(n), 'pfA'] = \
                                        gdf.loc[gdf['stream'] == visiting, 'pfA'].values.astype(int)[0]
                                        tovisit = tovisit.append(n, ignore_index=True)
                                        tovisit_temp = tovisit.to_list()
                                        tovisit_temp.pop(0)
                                        tovisit = pd.Series(tovisit_temp)
                                    else:  # main junction (main tributary gets +1 [even number], main channel gets +2 [odd number])
                                        gdf.loc[gdf['stream'].isin(n), 'pfA'] = \
                                        gdf.loc[gdf['stream'] == visiting, 'pfA'].values.astype(int)[0] + 2 * 10 ** (
                                                    iHack_max - j)
                                        gdf.loc[gdf['stream'].isin(a), 'pfA'] = \
                                        gdf.loc[gdf['stream'] == visiting, 'pfA'].values.astype(int)[0] + 1 * 10 ** (
                                                    iHack_max - j)
                                        tovisit = tovisit.append(n, ignore_index=True)
                                        tovisit_temp = tovisit.to_list()
                                        tovisit_temp.pop(0)
                                        tovisit = pd.Series(tovisit_temp)
                                else:  # source stream, go back to junction
                                    tovisit_temp = tovisit.to_list()
                                    tovisit_temp.pop(0)
                                    tovisit = pd.Series(tovisit_temp)
                                i = i + 1  # number of iterations

                        bb = gdf['pfA'].unique()
                        for kk in np.arange(0, bb.size, 1):  # loop over subbasins (to explore in the next step)
                            gdf_tmp = gdf.loc[gdf['pfA'] == bb[kk]]
                            gdf_tmp.loc[~gdf_tmp['next_strea'].isin(gdf_tmp['stream']), 'next_strea'] = -1  # outlet
                            for i in np.arange(1, gdf.shape[1] - iNum_col):
                                gdf_tmp.loc[~gdf_tmp['prev_str0' + str(i)].isin(gdf_tmp['stream']), 'prev_str0' + str(
                                    i)] = 0  # source stream
                            gdf.loc[gdf['pfA'] == bb[kk]] = gdf_tmp

                        gdf_final.loc[gdf_final['pfA'] == aa[k]] = gdf

            if y == 0:
                gdf_final['MBID_new'] = gdf_final['MBID']
            else:
                mbid_new = cc.max() + 1
                gdf_final['MBID_new'] = mbid_new
                cc = np.append(cc, mbid_new)

            gdf_final = gdf_final.dropna(subset=['pfA'])
            gdf_orig_mb.loc[gdf_orig_mb['stream'].isin(gdf_final['stream'])] = gdf_final

        gdf_orig.loc[gdf_orig['MBID'] == cc[m]] = gdf_orig_mb

    gdf_orig['pfA'] = gdf_orig['pfA'].astype(int)
    gdf_orig['MBID_new'] = gdf_orig['MBID_new'].astype(int)
    gdf_orig['FiscalCode'] = sDomain + '_' + gdf_orig['MBID_new'].astype(str) + '_' + gdf_orig['pfA'].astype(str)
    gdf_orig['next_strea'] = gdf_bkp[
        'next_strea']  # setting original next_streams and previous streams from the backup gdf
    for i in np.arange(1, gdf_orig.shape[1] - iNum_col):
        gdf_orig['prev_str0' + str(i)] = gdf_bkp['prev_str0' + str(i)]

    # export to shp
    sFileNameWeOutput = sInput[len(sPathInput) + 1:-5] + str('_pf_tmp')
    sFileNameOutput = '%s.shp' % (str(sFileNameWeOutput))
    sFileNameOutput_prj = '%s.prj' % (str(sFileNameWeOutput))
    sOutput = os.path.join(sPathInput, sFileNameOutput)
    schema = gpd.io.file.infer_schema(gdf_orig)
    gdf_orig.to_file(sOutput, driver='ESRI Shapefile')

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(input_epsg))
    srs.MorphToESRI()
    streams_mod_prjfile = os.path.join(sPathInput, sFileNameOutput_prj)
    file_prj = open(streams_mod_prjfile, 'w')
    file_prj.write(srs.ExportToWkt())
    file_prj.close()

    # Reorder columns in table of contents
    sFileNameWeOutput2 = sInput[len(sPathInput) + 1:-5] + str('_pf')
    sFileNameOutput2 = '%s.shp' % (str(sFileNameWeOutput2))
    sFileNameOutput_prj2 = '%s.prj' % (str(sFileNameWeOutput2))
    sOutput2 = os.path.join(sPathOutput, sFileNameOutput2)
    ogr2ogr_command = 'ogr2ogr %s %s -sql " SELECT MBID AS MBID, cat AS cat, cum_length AS cum_length, elev_drop AS elev_drop, flow_accum AS flow_accum, gradient AS gradient, hack AS hack, horton AS horton, length AS lenght, next_strea AS next_strea, out_dist AS out_dist, out_drop AS out_drop, outlet_ele AS outlet_ele, prev_str01 AS prev_str01, prev_str02 AS prev_str02, scheidegge AS scheidegge, shreve AS shreve, source_ele AS source_ele, strahler AS strahler, stream AS stream, cat_bas AS cat_bas, pfA AS pfA, FIscalCode AS FiscalCode, MBID_new AS MBID_new from %s"' % (
    sOutput2, sOutput, sFileNameWeOutput)
    print('-> Executing GDAL ogr2ogr..')
    print('%s' % ogr2ogr_command)
    os.system(ogr2ogr_command)
    #p = subprocess.Popen(ogr2ogr_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    #                     env={'GDAL_DATA': os.environ['GDAL_DATA']})
    #out, err = p.communicate()
    #if err:
    #    print('-> GDAL ogr2ogr Error: %s' % err)
    #else:
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(input_epsg))
    srs.MorphToESRI()
    streams_mod_prjfile = os.path.join(sPathOutput, sFileNameOutput_prj2)
    file_prj = open(streams_mod_prjfile, 'w')
    file_prj.write(srs.ExportToWkt())
    file_prj.close()

    return sOutput2

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to find all upstream river stretches
def Pfaf_find_upstream(pfA, maxpfA):
    lowlim = pfA
    uplim = pfA
    for i in range(len(str(pfA)) - 1, -1, -1):

        d = int(str(pfA)[i])
        if d == 0:
            continue
        else:
            if d % 2 == 0:
                # print('basin')
                if 'uplim' not in locals():
                    uplim = pfA
                break
            else:
                # print('interbasin')
                if i == 0:
                    uplim = maxpfA
                else:
                    uplim = int(str(int(str(pfA)[0:i]) + 1) + '0' * (len(str(pfA)) - i))

    return uplim

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to find a maximum of 2 upstream stretches
def Pfaf_find_Nextupstreams(pfaf_tmp, tmp_Upstream, righe, cat_tmp, head_basins, pd_DB):
    tmp_stream00 = tmp_Upstream.loc[(tmp_Upstream['pfA'] == pfaf_tmp)]

    prev_stream01 = tmp_Upstream.loc[(tmp_Upstream['pfA'] == pfaf_tmp)]['prev_str01']
    prev_stream02 = tmp_Upstream.loc[(tmp_Upstream['pfA'] == pfaf_tmp)]['prev_str02']

    if prev_stream01.values[0] > 0:
        tmp_stream01 = tmp_Upstream.loc[(tmp_Upstream['stream'] == prev_stream01.values[0])]
        # tmp_stream01_a = tmp_Upstream.loc[(tmp_Upstream['stream'] == tmp_stream01['prev_str01'].values[0])]
        # tmp_stream01_b = tmp_Upstream.loc[(tmp_Upstream['stream'] == tmp_stream01['prev_str02'].values[0])]
        if prev_stream02.values[0] > 0:
            tmp_stream02 = tmp_Upstream.loc[(tmp_Upstream['stream'] == prev_stream02.values[0])]
            # tmp_stream02_a = tmp_Upstream.loc[(tmp_Upstream['stream'] == tmp_stream02['prev_str01'].values[0])]
            # tmp_stream02_b = tmp_Upstream.loc[(tmp_Upstream['stream'] == tmp_stream02['prev_str02'].values[0])]
            tmp_Upstream_new = pd.concat([tmp_stream00, tmp_stream01, tmp_stream02],
                                         ignore_index=True)  # , tmp_stream01_a, tmp_stream01_b, tmp_stream02_a, tmp_stream02_b],
        # ignore_index=True)

        else:
            print('caso che non dovrebbe esistere!!  pfA={})'.format(pfaf_tmp))
            sys.exit()
    else:
        tmp_Upstream_new = tmp_Upstream.loc[(tmp_Upstream['pfA'] == pfaf_tmp)]  # ultima asta a monte
    # next_cat = pd_DB.loc[(pd_DB['stream'] == righe['next_strea'])]['cat_bas']
    # head_basins.append(cat_tmp)
    # head_basins.append(next_cat.values[0])

    return tmp_Upstream_new, head_basins

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to find a maximum of 2 downstream stretches
def Pfaf_find_downstream(pfA, minpfA):
    lowlim = pfA
    uplim = pfA
    for i in range(len(str(pfA)) - 1, -1, -1):

        d = int(str(pfA)[i])
        if d == 0:
            continue
        elif d == 1:
            continue
        else:
            if d % 2 == 0:
                if i == 0:
                    lowlim = int(str(int(str(pfA)[i]) - 1) + '0' * (len(str(pfA)) - i - 1))
                    uplim = pfA
                else:
                    lowlim = int(str(int(str(pfA)[0:i])) + str(int(str(pfA)[i]) - 1) + '0' * (len(str(pfA)) - i - 1))
                    uplim = pfA
                    break
            else:
                # print('interbasin')
                if i == 0:
                    lowlim = minpfA
                    # uplim = pfA
                    uplim = int(str(int(str(pfA)[i]) - 1) + '0' * (len(str(pfA)) - i - 1))

                else:
                    lowlim = int(str(int(str(pfA)[0:i])) + str(int(str(pfA)[i]) - 2) + '0' * (len(str(pfA)) - i - 1))
                    # uplim = pfA
                    uplim = int(str(int(str(pfA)[0:i])) + str(int(str(pfA)[i]) - 1) + '0' * (len(str(pfA)) - i - 1))
                    break

    return lowlim, uplim


# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to explore upstream junctions
def explore_tree(db_str, num_str, num_col):
    num_up = []
    for i in np.arange(1, db_str.shape[1] - num_col):
        num_up_tmp = db_str["prev_str0" + str(i)].loc[db_str["stream"] == num_str].squeeze()
        num_up.append(num_up_tmp)
    num_up_ser = pd.Series(num_up, name="stream")
    return num_up_ser

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to extract junction streams
def get_junction_streams(db_str):
    dimen = db_str.shape[0]
    pos_jun = []
    for k in np.arange(0, dimen, 1):
        num_dwn = db_str["next_strea"].to_list()[k]
        if num_dwn == -1:
            delta_hack = 1
        else:
            hack_cur = db_str["hack"].to_list()[k]
            if db_str["hack"][db_str["stream"] == num_dwn].empty == False:
                hack_dwn = db_str["hack"][db_str["stream"] == num_dwn].squeeze()
            else:
                hack_dwn = 0
            delta_hack = hack_cur - hack_dwn
        if delta_hack == 1:
            pos_jun.append(k)
    db_jun = db_str.iloc[pos_jun]
    return db_jun

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to check geometry
def isvalid(geom):  # Check Geometry
    try:
        shape(geom)
        return 1
    except:
        return 0

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to select 4 main tributaries
def select_mainstreams(db_str, hack_num):
    hack = db_str.loc[db_str["hack"] == hack_num]
    hack_ord = hack.sort_values("flow_accum", ascending=False)
    mainstreams = hack_ord["stream"][0:4]
    return mainstreams
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------