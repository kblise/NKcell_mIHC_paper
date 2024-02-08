# -*- coding: utf-8 -*-
"""
Code author: Katie E. Blise
Date: February 2024

This .py file contains all code needed to reproduce the spatial analyses and figures in the paper:
"Natural Killer cells occupy unique spatial neighborhoods in human HER2- and HER2+ breast cancers"

This .py file contains multiple functions:
    
    ***FUNCTIONS FOR FUNCTIONAL PROXIMITY ANALYSES***
    getCsvList() = list of all mIHC files for all analyses
    nkFunTumSpatial() = gets functional status of NK cells and spatial proximity to neoplastic tumor cells
    tumFunNKSpatial() = gets functional status of neoplastic tumor cells and spatial proximity to NK cells
    
    ***FUNCTIONS FOR NEIGHBORHOOD ANALYSES***
    makeNeighborhoods() = calculates spatial neighbors of seed cells within set distance
    elbowMethod() = runs elbow method to determine optimal number of clusters
    clusterNeighborhoods() = clusters neighborhoods based upon cellular compositions
    createCsvsWithClusterCol = creates new mIHC csvs with cluster column denoting NK cell neighborhood assignment 
    clusterCountPerROI = calculates how many seed cells are assigned to each cluster per patient, ROI 
   
    ***FUNCTIONS TO GENERATE RESULTS (which call to the above functions)***
    fig3() = generates Figures 3A-3D, Supplementary Figures S5A-S5C
    fig4() = generates Figures 4C-4E
    fig5() = generates Figures 5B-5F, Supplementary Figures S7A-S7C
   
Note: This program assumes the following items live in the same directory as this .py file:
    - 'data' folder, which houses 2 folders:
        -'mIHC_files' folder, which houses all mIHC data (.csv files)
        -'metadata' folder, which houses clinical data
    -'results' folder, which houses 2 folders:
        -'dfCreated' folder, which will store dataframes created by this code, and houses 1 folder:
            -'updatedCsvs' folder within 'dfCreated' folder, which will store revised mIHC csv files with neighborhood clustering assignments
        -'figures' folder, which will store figures created by this code

This program is intended for Python version 3.
"""



def getCsvList():
    '''
    This function generates a list of all mIHC files in specific order for reproducible results (when running neighborhood clustering analyses)
    Input parameters:
        None
    Outputs:
        returns: csvList = orgnaized list of all mIHC files in the dataset
    '''
    csvList = ['M10_TT1241A_ROI02', 'M10_TT1241A_ROI03', 'M10_TT1241A_ROI04', 'M10_TT1241A_ROI05', 'M10_TT1241A_ROI06', 
               'M10_TT1241A_ROI07', 'M10_TT1241A_ROI08', 'M10_TT1241A_ROI09', 'M10_TT1241A_ROI10', 'M10_TT1241A_ROI11', 
               'M11_TT1020A_ROI01', 'M11_TT1020A_ROI02', 'M11_TT1020A_ROI04', 'M11_TT1020A_ROI05', 'M11_TT1020A_ROI06', 
               'M11_TT1020A_ROI07', 'M11_TT1020A_ROI08', 'M11_TT1020A_ROI10', 'M11_TT1020A_ROI11', 'M11_TT1020A_ROI12', 
               'M11_TT1020A_ROI13', 'M12_TT1619A_ROI01', 'M12_TT1619A_ROI02', 'M12_TT1619A_ROI03', 'M12_TT1619A_ROI04', 
               'M12_TT1619A_ROI05', 'M12_TT1619A_ROI07', 'M12_TT1619A_ROI08', 'M13_TT1604A_ROI01', 'M13_TT1604A_ROI02', 
               'M13_TT1604A_ROI03', 'M13_TT1604A_ROI04', 'M13_TT1604A_ROI05', 'M13_TT1604A_ROI06', 'M13_TT1604A_ROI07', 
               'M14_TT1606A_ROI01', 'M14_TT1606A_ROI02', 'M14_TT1606A_ROI03', 'M14_TT1606A_ROI04', 'M14_TT1606A_ROI05', 
               'M14_TT1606A_ROI06', 'M14_TT1606A_ROI07', 'M14_TT1606A_ROI08', 'M14_TT1606A_ROI09', 'M14_TT1606A_ROI10', 
               'M14_TT1606A_ROI11', 'M15_TT1377A_ROI01', 'M15_TT1377A_ROI02', 'M15_TT1377A_ROI03', 'M15_TT1377A_ROI04', 
               'M15_TT1377A_ROI05', 'M15_TT1377A_ROI06', 'M15_TT1377A_ROI07', 'M15_TT1377A_ROI08', 'M15_TT1377A_ROI09', 
               'M15_TT1377A_ROI10', 'M15_TT1377A_ROI11', 'M16_TT1109A_ROI01', 'M16_TT1109A_ROI02', 'M16_TT1109A_ROI04', 
               'M16_TT1109A_ROI05', 'M16_TT1109A_ROI06', 'M16_TT1109A_ROI08', 'M16_TT1109A_ROI09', 'M16_TT1109A_ROI10', 
               'M17_TT1408A_ROI01', 'M17_TT1408A_ROI02', 'M17_TT1408A_ROI03', 'M17_TT1408A_ROI04', 'M17_TT1408A_ROI05', 
               'M17_TT1408A_ROI06', 'M17_TT1408A_ROI07', 'M17_TT1408A_ROI08', 'M17_TT1408A_ROI09', 'M17_TT1408A_ROI10', 
               'M18_TT1320A_ROI01', 'M18_TT1320A_ROI02', 'M18_TT1320A_ROI03', 'M18_TT1320A_ROI04', 'M18_TT1320A_ROI05', 
               'M18_TT1320A_ROI06', 'M18_TT1320A_ROI07', 'M18_TT1320A_ROI08', 'M18_TT1320A_ROI09', 'M18_TT1320A_ROI10', 
               'M19_TT1249A_ROI01', 'M19_TT1249A_ROI02', 'M19_TT1249A_ROI03', 'M19_TT1249A_ROI04', 'M19_TT1249A_ROI05', 
               'M19_TT1249A_ROI06', 'M19_TT1249A_ROI07', 'M19_TT1249A_ROI08', 'M19_TT1249A_ROI09', 'M19_TT1249A_ROI10', 
               'M1_TT1305A_ROI01', 'M1_TT1305A_ROI02', 'M1_TT1305A_ROI03', 'M1_TT1305A_ROI04', 'M1_TT1305A_ROI05', 
               'M1_TT1305A_ROI06', 'M1_TT1305A_ROI07', 'M1_TT1305A_ROI08', 'M20_TT1104A_ROI01', 'M20_TT1104A_ROI02', 
               'M20_TT1104A_ROI03', 'M20_TT1104A_ROI04', 'M20_TT1104A_ROI05', 'M20_TT1104A_ROI06', 'M20_TT1104A_ROI07', 
               'M20_TT1104A_ROI08', 'M20_TT1104A_ROI09', 'M21_TT1409A_ROI01', 'M21_TT1409A_ROI02', 'M21_TT1409A_ROI03', 
               'M21_TT1409A_ROI04', 'M21_TT1409A_ROI05', 'M21_TT1409A_ROI06', 'M21_TT1409A_ROI08', 'M21_TT1409A_ROI09', 
               'M21_TT1409A_ROI10', 'M21_TT1409A_ROI11', 'M22_TT1303A_ROI01', 'M22_TT1303A_ROI02', 'M22_TT1303A_ROI03', 
               'M22_TT1303A_ROI04', 'M22_TT1303A_ROI05', 'M22_TT1303A_ROI06', 'M22_TT1303A_ROI07', 'M22_TT1303A_ROI08', 
               'M23_TT1202A_ROI01', 'M23_TT1202A_ROI02', 'M23_TT1202A_ROI03', 'M23_TT1202A_ROI04', 'M23_TT1202A_ROI05', 
               'M23_TT1202A_ROI06', 'M23_TT1202A_ROI07', 'M23_TT1202A_ROI08', 'M23_TT1202A_ROI09', 'M23_TT1202A_ROI10', 
               'M24_TT1113A_ROI02', 'M24_TT1113A_ROI03', 'M24_TT1113A_ROI04', 'M24_TT1113A_ROI05', 'M24_TT1113A_ROI06', 
               'M24_TT1113A_ROI07', 'M24_TT1113A_ROI08', 'M24_TT1113A_ROI09', 'M24_TT1113A_ROI10', 'M24_TT1113A_ROI11', 
               'M24_TT1113A_ROI12', 'M25_TT1205A_ROI01', 'M25_TT1205A_ROI02', 'M25_TT1205A_ROI03', 'M25_TT1205A_ROI04', 
               'M25_TT1205A_ROI05', 'M25_TT1205A_ROI06', 'M25_TT1205A_ROI07', 'M25_TT1205A_ROI08', 'M25_TT1205A_ROI09', 
               'M25_TT1205A_ROI10', 'M25_TT1205A_ROI11', 'M26_TT1616A_ROI01', 'M26_TT1616A_ROI02', 'M26_TT1616A_ROI03', 
               'M26_TT1616A_ROI04', 'M26_TT1616A_ROI05', 'M26_TT1616A_ROI06', 'M26_TT1616A_ROI07', 'M26_TT1616A_ROI08', 
               'M27_TT1120A_ROI01', 'M27_TT1120A_ROI02', 'M27_TT1120A_ROI03', 'M27_TT1120A_ROI04', 'M27_TT1120A_ROI05', 
               'M27_TT1120A_ROI06', 'M27_TT1120A_ROI07', 'M28_TT1609A_ROI02', 'M28_TT1609A_ROI03', 'M28_TT1609A_ROI04', 
               'M28_TT1609A_ROI05', 'M28_TT1609A_ROI06', 'M28_TT1609A_ROI07', 'M28_TT1609A_ROI08', 'M28_TT1609A_ROI09', 
               'M28_TT1609A_ROI10', 'M28_TT1609A_ROI11', 'M28_TT1609A_ROI12', 'M29_TT1217A_ROI01', 'M29_TT1217A_ROI02', 
               'M29_TT1217A_ROI03', 'M29_TT1217A_ROI04', 'M29_TT1217A_ROI05', 'M29_TT1217A_ROI06', 'M29_TT1217A_ROI07', 
               'M29_TT1217A_ROI09', 'M29_TT1217A_ROI10', 'M2_TT1316A_ROI01', 'M2_TT1316A_ROI02', 'M2_TT1316A_ROI03', 
               'M2_TT1316A_ROI04', 'M2_TT1316A_ROI05', 'M2_TT1316A_ROI06', 'M2_TT1316A_ROI07', 'M2_TT1316A_ROI08', 
               'M2_TT1316A_ROI09', 'M30_TT1417A_ROI01', 'M30_TT1417A_ROI02', 'M30_TT1417A_ROI03', 'M30_TT1417A_ROI04', 
               'M30_TT1417A_ROI05', 'M30_TT1417A_ROI06', 'M30_TT1417A_ROI07', 'M30_TT1417A_ROI08', 'M30_TT1417A_ROI09', 
               'M30_TT1417A_ROI10', 'M30_TT1417A_ROI11', 'M3_TT1114A_ROI01', 'M3_TT1114A_ROI02', 'M3_TT1114A_ROI03', 
               'M3_TT1114A_ROI04', 'M3_TT1114A_ROI05', 'M3_TT1114A_ROI06', 'M3_TT1114A_ROI07', 'M3_TT1114A_ROI08', 
               'M4_TT1306A_ROI01', 'M4_TT1306A_ROI02', 'M4_TT1306A_ROI03', 'M4_TT1306A_ROI04', 'M4_TT1306A_ROI05', 
               'M4_TT1306A_ROI06', 'M4_TT1306A_ROI07', 'M4_TT1306A_ROI08', 'M5_TT1009A_ROI01', 'M5_TT1009A_ROI02', 
               'M5_TT1009A_ROI04', 'M6_TT1206A_ROI01', 'M6_TT1206A_ROI02', 'M6_TT1206A_ROI03', 'M6_TT1206A_ROI04', 
               'M6_TT1206A_ROI05', 'M6_TT1206A_ROI06', 'M6_TT1206A_ROI07', 'M6_TT1206A_ROI08', 'M6_TT1206A_ROI09', 
               'M6_TT1206A_ROI10', 'M7_TT1322A_ROI01', 'M7_TT1322A_ROI02', 'M7_TT1322A_ROI03', 'M7_TT1322A_ROI04', 
               'M7_TT1322A_ROI05', 'M7_TT1322A_ROI06', 'M7_TT1322A_ROI07', 'M8_TT1502A_ROI02', 'M8_TT1502A_ROI03', 
               'M8_TT1502A_ROI04', 'M8_TT1502A_ROI06', 'M8_TT1502A_ROI07', 'M8_TT1502A_ROI08', 'M8_TT1502A_ROI09', 
               'M9_TT1412A_ROI01', 'M9_TT1412A_ROI02', 'M9_TT1412A_ROI03', 'M9_TT1412A_ROI04', 'M9_TT1412A_ROI05', 
               'M9_TT1412A_ROI06', 'M9_TT1412A_ROI07', 'M9_TT1412A_ROI08', 'D12_BB1944A_ROI01', 'D12_BB1944A_ROI02', 
               'D13_TN1248A_ROI01', 'D14_BB2043A_ROI01', 'D14_BB2043A_ROI02', 'D14_BB2043A_ROI03', 'D15_BB2069A_ROI01', 
               'D15_BB2069A_ROI02', 'D15_BB2069A_ROI03', 'D15_BB2069A_ROI04', 'D16_BB2014A_ROI01', 'D16_BB2014A_ROI02', 
               'D16_BB2014A_ROI03', 'D16_BB2014A_ROI04', 'D17_TN1253A_ROI01', 'D17_TN1253A_ROI02', 'D19_BB2090A_ROI02', 
               'D19_BB2090A_ROI03', 'D20_BB2044A_ROI01', 'D20_BB2044A_ROI02', 'D20_BB2044A_ROI03', 'D21_BR1348A_ROI01', 
               'D21_BR1348A_ROI02', 'D24_BR1345A_ROI01', 'D24_BR1345A_ROI02', 'D25_TN1220A_ROI01', 'D25_TN1220A_ROI02', 
               'D26_TN1242A_ROI01', 'D2_TN1233A_ROI01', 'D2_TN1233A_ROI02', 'D3_TN1247A_ROI01', 'D4_BB2028A_ROI01', 
               'D5_TN1252A_ROI01', 'D6_BB2088A_ROI01', 'D6_BB2088A_ROI02', 'D6_BB2088A_ROI03', 'D6_BB2088A_ROI04', 
               'D6_BB2088A_ROI05', 'D6_BB2088A_ROI06', 'D7_TN1246A_ROI01', 'D7_TN1246A_ROI02', 'D7_TN1246A_ROI03', 
               'D10_BB2006A_ROI01', 'D10_BB2006A_ROI02', 'D11_BB2120A_ROI01', 'D11_BB2120A_ROI02', 'D11_BB2120A_ROI03', 
               'D11_BB2120A_ROI04', 'D18_BB2096A_ROI01', 'D1_TN1255A_ROI01', 'D22_BB2091A_ROI01', 'D22_BB2091A_ROI02', 
               'D22_BB2091A_ROI03', 'D23_BB2172A_ROI01', 'D23_BB2172A_ROI02', 'D23_BB2172A_ROI03', 'D23_BB2172A_ROI04', 
               'D23_BB2172A_ROI05', 'D23_BB2172A_ROI06', 'D8_BB2127A_ROI01', 'D8_BB2127A_ROI02', 'D8_BB2127A_ROI03', 
               'D9_BB2095A_ROI01', 'D9_BB2095A_ROI02', 'D9_BB2095A_ROI03', 'D9_BB2095A_ROI04', 'D9_BB2095A_ROI05']

    return csvList



def nkFunTumSpatial(path,csvList,distThresh):
    '''
    This function identifies the functional status of NK cells that are proximal and distal to neoplastic epithelial cells
    Input parameters:
        path = cwd
        csvList = list of mIHC files in the dataset
        distThresh = distance to stratify proximal vs distal (in px); note 1 µm = 2 px
    Outputs:
        Saves one csv with proportion of NK cells expressing functional markers that are proximal vs distal to neoplastic cells per patient. Csv saved to the /results/dfCreated/ folder.
    '''
    
    import pandas as pd
    from scipy import spatial

    cellsToKeep = ['Tumor cells','CD56+ NKP46+ NK','CD56+ NKP46- NK','CD56- NKP46+ NK']
    neighType = 'Tumor cells'
    seedList = ['CD56+ NKP46+ NK','CD56+ NKP46- NK','CD56- NKP46+ NK']

    #empty lists to store functional results - for NK seeds with a tumor neighbor
    fileList = []
    seedIdxList = []
    cd16List = []
    cd57List = []
    ki67List = []
    nkg2dList = []
    pd1List = []
    tim3List = []
    grzbList = []
    
    #empty lists to store functional results - for NK seeds withOUT a tumor neighbor
    fileFarList = []
    seedIdxFarList = []
    cd16FarList = []
    cd57FarList = []
    ki67FarList = []
    nkg2dFarList = []
    pd1FarList = []
    tim3FarList = []
    grzbFarList = []
    
    #loop through each file in the csvList
    for file in csvList:
        #read original csv
        df = pd.read_csv(path+'/data/mIHC_files/'+file+'.csv', index_col=0)

        #only care about specific cells so filter down to speed up neighbor loop
        filt_df = df[(df['class'].isin(cellsToKeep))]

        ##get nearest neighbors of seed cells defined by seed param
        #create np array of just x,y coordinates
        ptsArray = filt_df[['Location_Center_X','Location_Center_Y']].values

        #create kdtree
        tree = spatial.KDTree(ptsArray)

        #loop through each cell in the filt_df and check its neighbors if it's an NK cell
        for i in range(len(ptsArray)):

            classType = filt_df['class'].values[i]

            #only check neighbors if the seed cell is the desired classType aka an NK cell
            if classType in seedList:
                neighs = tree.query_ball_point(ptsArray[i], distThresh) #get neighbors within distance

                #check if there is a tumor cell as a neighbor within the distance; loop through each neighbor and get its class
                #start as not having a tumor cell as a neighbor
                tumNeigh = False
                while tumNeigh == False:
                    for j in neighs:

                        #get its class
                        neighClass = filt_df['class'].values[j]

                        #if it's a tumor cell, get fun status of seed, eventually change to True
                        if neighClass == neighType:

                            #get seed cell index (original df.loc index), function
                            seedIdx = filt_df.iloc[i].name
                            cd16 = filt_df['Cellsp_CD16p'].values[i]
                            cd57 = filt_df['Cellsp_CD57p'].values[i]
                            ki67 = filt_df['Cellsp_Ki67p'].values[i]
                            nkg2d = filt_df['Cellsp_NKG2Dp'].values[i]
                            pd1 = filt_df['Cellsp_PD1p'].values[i]
                            tim3 = filt_df['Cellsp_TIM3p'].values[i]
                            grzb = filt_df['Cellsp_GRZBp'].values[i]

                            #add values to lists to store in df later
                            fileList.append(file)
                            seedIdxList.append(seedIdx)
                            cd16List.append(cd16)
                            cd57List.append(cd57)
                            ki67List.append(ki67)
                            nkg2dList.append(nkg2d)
                            pd1List.append(pd1)
                            tim3List.append(tim3)
                            grzbList.append(grzb)
                            
                            tumNeigh = True
                            break

                        #if not a Tumor cell neighbor, continue checking the next neighbor
                        else:
                            continue

                    #if tumNeigh still false after looping through all neighbors, the seed cell is not close to tumor cells
                    #get its function too - store in far df
                    if tumNeigh == False:
                        seedIdxFar = filt_df.iloc[i].name
                        cd16Far = filt_df['Cellsp_CD16p'].values[i]
                        cd57Far = filt_df['Cellsp_CD57p'].values[i]
                        ki67Far = filt_df['Cellsp_Ki67p'].values[i]
                        nkg2dFar = filt_df['Cellsp_NKG2Dp'].values[i]
                        pd1Far = filt_df['Cellsp_PD1p'].values[i]
                        tim3Far = filt_df['Cellsp_TIM3p'].values[i]
                        grzbFar = filt_df['Cellsp_GRZBp'].values[i]
                        
                        #add values to lists to store in df later
                        fileFarList.append(file)
                        seedIdxFarList.append(seedIdxFar)
                        cd16FarList.append(cd16Far)
                        cd57FarList.append(cd57Far)
                        ki67FarList.append(ki67Far)
                        nkg2dFarList.append(nkg2dFar)
                        pd1FarList.append(pd1Far)
                        tim3FarList.append(tim3Far)
                        grzbFarList.append(grzbFar)

                        #then make tumNeigh True to exit while loop and move to the next seed       
                        tumNeigh = True

    #store results in df - for seeds with a Tumor neighbor
    dfFunClose = pd.DataFrame()
    dfFunClose['file'] = fileList
    dfFunClose['seedIdx'] = seedIdxList
    dfFunClose['cd16'] = cd16List
    dfFunClose['cd57'] = cd57List
    dfFunClose['ki67'] = ki67List
    dfFunClose['nkg2d'] = nkg2dList
    dfFunClose['pd1'] = pd1List
    dfFunClose['tim3'] = tim3List
    dfFunClose['grzb'] = grzbList
    
    #store results in df - for seeds withOUT a Tumor neighbor
    dfFunFar = pd.DataFrame()
    dfFunFar['file'] = fileFarList
    dfFunFar['seedIdx'] = seedIdxFarList
    dfFunFar['cd16'] = cd16FarList
    dfFunFar['cd57'] = cd57FarList
    dfFunFar['ki67'] = ki67FarList
    dfFunFar['nkg2d'] = nkg2dFarList
    dfFunFar['pd1'] = pd1FarList
    dfFunFar['tim3'] = tim3FarList
    dfFunFar['grzb'] = grzbFarList
 
    #now merge close and far dfs into one df and organize per patient
    locat = ['close','far']
    nameDict = {'close':dfFunClose,'far':dfFunFar}

    #empty lists to store % cells positive per patient
    cd16List = []
    cd57List = []
    ki67List = []
    nkg2dList = []
    pd1List = []
    tim3List = []
    grzbList = []
    totalList = [] #stores raw count of NK cells for that patient at that location
    her2List = []
    locatList = []
    ptList = []

    #read clinical df for HER2 status
    dfClin = pd.read_csv(path+'/data/metadata/clinicalData.csv',index_col=0)

    #loop through close then far dfs created earlier
    for l in locat:

        dfName = nameDict[l]

        #get just patient name
        dfName['file'] = dfName['file'].str.slice(0,-6)

        #subset to just one patient at a time
        for pt in dfName['file'].unique():

            dfPt = dfName[dfName['file'] == pt]

            #get total count of nk cells
            totNK = len(dfPt)

            #sum the functional markers for the patient - across all ROIs
            dfTotal = dfPt.sum(axis=0)

            #get patient's her2 status
            her2 = dfClin.loc[pt,'HER2'] #0 = HER2-, 1 = HER2+

            cd16List.append(dfTotal['cd16']/totNK*100)
            cd57List.append(dfTotal['cd57']/totNK*100)
            ki67List.append(dfTotal['ki67']/totNK*100)
            nkg2dList.append(dfTotal['nkg2d']/totNK*100)
            pd1List.append(dfTotal['pd1']/totNK*100)
            tim3List.append(dfTotal['tim3']/totNK*100)
            grzbList.append(dfTotal['grzb']/totNK*100)
            totalList.append(totNK)
            her2List.append(her2)
            locatList.append(l)
            ptList.append(pt)

    #put results into a df to plot
    dfFun = pd.DataFrame([cd16List,cd57List,ki67List,nkg2dList,pd1List,tim3List,grzbList,totalList,her2List,locatList,ptList]).T
    dfFun.columns = ['CD16','CD57','KI67','NKG2D','PD1','TIM3','GRZB','Total NK Cells','HER2','Location','Patient']

    #save dfFun to csv - this gets used to create figures
    dfFun.to_csv(path+'/results/dfCreated/dfNKFun_TumorSpatial_all'+str(distThresh)+'.csv')    
    
    
    
def tumorFunNKspatial(path,csvList,distThresh):
    '''
    This function identifies the functional status of neoplastic cells that are proximal and distal to NK cells
    Input parameters:
        path = cwd
        csvList = list of mIHC files in the dataset
        distThresh = distance to stratify proximal vs distal (in px); note 1 µm = 2 px
    Outputs:
        Saves one csv with proportion of neoplastic cells expressing functional markers that are proximal vs distal to NK cells per patient. Csv saved to the /results/dfCreated/ folder.
    '''

    import pandas as pd
    from scipy import spatial

    cellsToKeep = ['Tumor cells','CD56+ NKP46+ NK','CD56+ NKP46- NK','CD56- NKP46+ NK']
    seed = 'Tumor cells'
    nkList = ['CD56+ NKP46+ NK','CD56+ NKP46- NK','CD56- NKP46+ NK']

    #empty lists to store functional results - for seeds with a NK neighbor
    fileList = []
    seedIdxList = []
    hla1List = []
    ki67List = []
    pdl1List = []
    caixList = []

    #empty lists to store functional results - for seeds withOUT a NK neighbor
    fileFarList = []
    seedIdxFarList = []
    hla1FarList = []
    ki67FarList = []
    pdl1FarList = []
    caixFarList = []

    #loop through each file in the csvList
    for file in csvList:
        #read original csv
        df = pd.read_csv(path+'/data/mIHC_files/'+file+'.csv', index_col=0)

        #only care about specific cells so filter down to speed up neighbor loop
        filt_df = df[(df['class'].isin(cellsToKeep))]

        ##get nearest neighbors of seed cells defined by seed param
        #create np array of just x,y coordinates
        ptsArray = filt_df[['Location_Center_X','Location_Center_Y']].values

        #create kdtree
        tree = spatial.KDTree(ptsArray)

        #loop through each cell in the filt_df and check its neighbors if it's a tumor cell
        for i in range(len(ptsArray)):
            classType = filt_df['class'].values[i]

            #only check neighbors if the seed cell is the desired classType
            if classType == seed:

                neighs = tree.query_ball_point(ptsArray[i], distThresh) #get neighbors within distance
                #check if there is an NK cell as a neighbor within the distance; loop through each neighbor and get its class

                #start as not having an NK cell as a neighbor
                nkNeigh = False

                while nkNeigh == False:
                    for j in neighs:
                        #get its class
                        neighClass = filt_df['class'].values[j]

                        #if it's an NK cell, get fun status, eventually change to True
                        if neighClass in nkList:

                            #get seed cell index (original df.loc index), function
                            seedIdx = filt_df.iloc[i].name
                            hla1 = filt_df['Cellsp_HLAIIp'].values[i] #in csv it's listed as HLAII even though it's HLAI
                            ki67 = filt_df['Cellsp_Ki67p'].values[i]
                            pdl1 = filt_df['Cellsp_PDL1p'].values[i]
                            caix = filt_df['Cellsp_CAIXp'].values[i]

                            #add values to lists to store in df later
                            fileList.append(file)
                            seedIdxList.append(seedIdx)
                            hla1List.append(hla1)
                            ki67List.append(ki67)
                            pdl1List.append(pdl1)
                            caixList.append(caix)
                            nkNeigh = True
                            break

                        #if not a NK cell neighbor, continue checking the next neighbor
                        else:
                            continue

                    #if nkNeigh still false after looping through all neighbors, the seed cell is not close to NK cells
                    #get its function too - store in far df
                    if nkNeigh == False:
                        seedIdxFar = filt_df.iloc[i].name
                        hla1Far = filt_df['Cellsp_HLAIIp'].values[i] #in csv it's listed as HLAII even though it's HLAI
                        ki67Far = filt_df['Cellsp_Ki67p'].values[i]
                        pdl1Far = filt_df['Cellsp_PDL1p'].values[i]
                        caixFar = filt_df['Cellsp_CAIXp'].values[i]

                        #add values to lists to store in df later
                        fileFarList.append(file)
                        seedIdxFarList.append(seedIdxFar)
                        hla1FarList.append(hla1Far)
                        ki67FarList.append(ki67Far)
                        pdl1FarList.append(pdl1Far)
                        caixFarList.append(caixFar)

                        #then make nkNeigh True to exit while loop and move to the next seed       
                        nkNeigh = True

    #store results in df - for seeds with a NK neighbor
    dfFunClose = pd.DataFrame()
    dfFunClose['file'] = fileList
    dfFunClose['seedIdx'] = seedIdxList
    dfFunClose['hla1'] = hla1List
    dfFunClose['ki67'] = ki67List
    dfFunClose['pdl1'] = pdl1List
    dfFunClose['caix'] = caixList

    #store results in df - for seeds withOUT a NK neighbor
    dfFunFar = pd.DataFrame()
    dfFunFar['file'] = fileFarList
    dfFunFar['seedIdx'] = seedIdxFarList
    dfFunFar['hla1'] = hla1FarList
    dfFunFar['ki67'] = ki67FarList
    dfFunFar['pdl1'] = pdl1FarList
    dfFunFar['caix'] = caixFarList

    locat = ['close','far']
    nameDict = {'close':dfFunClose,'far':dfFunFar}
    
    #empty lists to store % cells positive per patient
    hla1List = []
    ki67List = []
    pdl1List = []
    caixList = []
    totalList = [] #stores raw count of NK cells for that patient at that location
    her2List = []
    locatList = []
    ptList = []
    
    #read clinical df for HER2 status
    dfClin = pd.read_csv(path+'/data/metadata/clinicalData.csv',index_col=0)
    
    #loop through close then far tumor cell csvs
    for l in locat:
        dfName = nameDict[l]
        
        #get just patient name
        dfName['file'] = dfName['file'].str.slice(0,-6)
    
        #subset to just one patient at a time
        for pt in dfName['file'].unique():
    
            dfPt = dfName[dfName['file'] == pt]
    
            #get total count of tumor cells
            totTum = len(dfPt)
    
            #sum the functional markers for the patient - across all ROIs
            dfTotal = dfPt.sum(axis=0)
    
            #get patient's her2 status
            her2 = dfClin.loc[pt,'HER2'] #0 = HER2-, 1 = HER2+
    
            hla1List.append(dfTotal['hla1']/totTum*100)
            ki67List.append(dfTotal['ki67']/totTum*100)
            pdl1List.append(dfTotal['pdl1']/totTum*100)
            caixList.append(dfTotal['caix']/totTum*100)
    
            totalList.append(totTum)
            her2List.append(her2)
            locatList.append(l)
            ptList.append(pt)
                
    #put results into a df to plot
    dfFun = pd.DataFrame([hla1List,ki67List,pdl1List,caixList,totalList,her2List,locatList,ptList]).T
    dfFun.columns = ['HLA1','KI67','PDL1','CAIX','Total Tumor Cells','HER2','Location','Patient']
    
    #save dfFun to csv
    dfFun.to_csv(path+'/results/dfCreated/dfTumorFun_NKspatial_all'+str(distThresh)+'.csv')
    
    
    

def makeNeighborhoods(path,csvList,seedList,distThresh):
    '''
    This function generates spatial neighborhoods for NK cells within a specified radius.
    Input parameters:
        path = cwd
        csvList = list of mIHC files in the dataset
        seedList = phenotypes to generate neighborhoods for
        distThresh = distance to set radius for spatial neighborhoods, in px, 2 px = 1 µm
    Outputs:
        saves one csv to 'dfCreated' folder with neighbors of each seed cell
    '''
        
    import pandas as pd
    from scipy import spatial

    #empty list to hold dictionaries to create new rows of dfClust; outside of for file in csvList loop
    allNeighList = []

    #loop through each ROI in csvList
    for file in csvList:

        #read df according to path
        df = pd.read_csv(path+'/data/mIHC_files/'+file+'.csv', index_col=0)

        #create filtered dataframe without noise or 'other cells'
        filt_df = df[(df['class'] != 'Other cells') & (df['class'] != 'Noise') ]

        #get all possible class values; for later counting - from prior knowledge of possible cell types
        classOptions = ['CD11B+ DCs',
                 'CD11B- CD68+ cells',
                 'CD11B- DCs',
                 'CD4 T cells',
                 'CD56+ NKP46+ NK',
                 'CD56+ NKP46- NK',
                 'CD56- NKP46+ NK',
                 'CD8 T cells',
                 'Myeloid other',
                 'Myelomonocytic cells',
                 'Other CD45+ cells',
                 'Tumor cells']
        
        #generate count variable names for final created df's columns
        classCountNames = []
        for c in classOptions: #loop through all possible neighboring cells
            classCountNames.append('count'+c)

        ##get nearest neighbors of seed cells defined by seed param
        #create np array of just x,y coordinates
        ptsArray = filt_df[['Location_Center_X','Location_Center_Y']].values

        #create kdtree
        tree = spatial.KDTree(ptsArray)

        #loop through each cell and check its neighbors if it's the right seed
        for i in range(len(ptsArray)):

            classType = filt_df['class'].values[i]

            #only check neighbors if the seed cell is the desired classType
            if classType in seedList:
                neighs = tree.query_ball_point(ptsArray[i], distThresh) #get neighbors within distance

                if i in neighs:
                    neighs.remove(i) #don't include itself as a neighbor

                neighClassList = [] #empty list to store neighboring cells' classes

                #loop through each neighbor and get its class; add it to neighClassList
                for j in neighs:
                    #get its class
                    neighClass = filt_df['class'].values[j]
                    neighClassList.append(neighClass)

                #get counts of neighboring cell types
                classCounts = []
                for c in classOptions: #loop through all possible neighboring cells
                    count = neighClassList.count(c) #count the specified cell type
                    classCounts.append(count) #add the counts of the specified cell type to classCounts

                #reset dictionary to hold counts of neighbors; one per seed
                seedDict = {}

                #add counts to a dictionary (one per seed cell); also add original seed cell's idx value 
                seedDict['file'] = file
                seedDict['index'] = filt_df.iloc[i].name #original seed cell's idx value in the form of its original csv's name aka df.loc index

                #for each possible neighboring cell class, add its count to seedDict both raw and as a %
                for n in range(len(classCountNames)):
                    seedDict[classCountNames[n]] = classCounts[n] #raw count

                    if sum(classCounts) != 0:
                        seedDict[classCountNames[n]+'%'] = classCounts[n]/sum(classCounts) #percentage
                    else: #avoid division by zero if there are no neighbors
                        seedDict[classCountNames[n]+'%'] = 0 #set % to zero if there are no neighbors

                #add each seed's neighbor dictionary to the overall list; one dictionary per row of df
                allNeighList.append(seedDict)

    #create one new df to hold data for clustering; pass in allNeighList as the data; format is one row per seed cell across all csvs
    #column names from a seedDict's keys (all seedDicts have the same keys)
    dfClust = pd.DataFrame(data = allNeighList, columns = list(seedDict.keys()))
    #NOTE: only using the columns (cells present) from the final csv that was analyzed. If a cell is not present at all in this csv then its column will not be present in the final csv created.

    #convert any NaN values to zeros [note that NaN values arise when a csv lacks any of a cell type that does exist in other csvs]
    dfClust = dfClust.fillna(0)

    #store dfClust as a csv
    dfClust.to_csv(path+'/results/dfCreated/dfNeighborhoodClusterNK'+str(distThresh)+'.csv')    
    
    

def elbowMethod(path,file,steps,save):

    '''
    This function runs the Elbow Method to determine the optimal number of clusters for k-means clustering.
    Input parameters:
        path = cwd
        file = name of file to run clustering on 
        steps = max number of clusters (k) to test
        save = if True, creates and saves the plot
        
    Output:
        optionally saves plot to 'figures' folder
    '''

    import pandas as pd
    from matplotlib import pyplot as plt
    from sklearn.cluster import MiniBatchKMeans #minibatchkmeans is better when n > 10,000 samples


    df = pd.read_csv(path+'/results/dfCreated/'+file+'.csv', index_col=0)
    #drop all rows that have no cells in the neighborhood (aka when the sum of count columns is zero)
    df['sum'] = df.iloc[:,2:].sum(axis=1)

    df = df[df['sum'] != 0]

    #generate column list to cluster on based on if there is a % in the column name
    colList = list(df.columns[['%' in col for col in list(df.columns)]])

    #get only features we want to cluster on
    df = df[colList]
    data = df.values

    #empty list to store error value
    wcss = []

    #calculate error for each k value (k=number of clusters)
    for k in range(1, steps):
        #generate kmeans model
        kmeans = MiniBatchKMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
        #fit model to data
        kmeans.fit(data)
        #add the sum of squares to wcss list; for plotting elbow
        wcss.append(kmeans.inertia_)

    #generate elbow plot and save (not results are not shown in manuscript)
    if save == True:
        plt.plot(range(1, steps), wcss)
        plt.title('Elbow Method')
        plt.xlabel('Number of clusters')
        plt.ylabel('WCSS')
        plt.savefig(path+'/results/figures/figure5_elbow_plot.png',format='png')
        plt.close()
    
   
def clusterNeighborhoods(path,file,k):
    '''
    This function runs k-means clustering on a given neighborhood clustering csv.
    The results are saved to a new csv.
    Input parameters:
        path = cwd
        file = name of file to run clustering on excluding the .csv
        k = number of clusters; use elbow method to determine optimal number
    Outputs:
        One csv is saved to the 'dfCreated/' folder with NK neighborhood cluster assignments
    '''

    import pandas as pd
    from sklearn.cluster import MiniBatchKMeans
    
    #read csv with neighborhood data    
    df = pd.read_csv(path+'/results/dfCreated/'+file+'.csv', index_col=0)

    #drop all rows that have no cells in the neighborhood (aka when the sum of count columns is zero)
    df['sum'] = df.iloc[:,2:].sum(axis=1)
    dfFilt =df[df['sum'] != 0]

    #get lists from original df to add back later after clustering
    roiList = list(dfFilt['file']) #get list of ROIs in order to add to dfNoNa later
    idxList = list(dfFilt['index']) #get list of cell indices in order to add to dfNoNa later

    #generate column list to cluster on based on if there is a % in the column name; this excludes the ki67, pdl1 columns too
    colList = list(dfFilt.columns[['%' in col for col in list(dfFilt.columns)]])

    #get only features we want to cluster on
    dfFilt = dfFilt[colList]
    data = dfFilt.values

    #=k-means clustering of cells with k clusters
    kmeans = MiniBatchKMeans(n_clusters=k, init='k-means++', max_iter=300, n_init=10, random_state=0)
    predict = kmeans.fit_predict(data) #fit model to data and predict index (cluster labels); same results as first fitting and then predicting

    #add predicted cluster labels to df as a new column and see how many cells are in each cluster
    dfFilt['cluster'] = predict
    dfFilt['cluster'].value_counts() #unecessary unless you want to see how many cells are in each cluster and then you should print it

    #add original ROI to df to check which ROIs are in each cluster
    dfFilt['file'] = roiList

    #need to add original indices to this df to check which ROIs are in each cluster; will also need to use ROI ID to pair with cell index (same index could be had by two cells from diff ROIs)
    dfFilt['index'] = idxList #idxList stores the row value of filt_df.iloc[row,column] command

    #save df to a csv
    dfFilt.to_csv(path+'/results/dfCreated/dfNeighClustered'+file[21:]+'k'+str(k)+'.csv')



def createCsvsWithClusterCol(path,name):
    '''
    This function takes a csv of all seed cells with their cluster IDs and outputs separate csvs for two hard-coded ROIs: D16_BB2014A_ROI01 and M27_TT1120A_ROI02
    The new ROIs contain a 'cluster' column which is populated with the corresponding cluster ID for each NK cell.
    Non-NK cells have their original cell classification ID populated in this column.

    Input parameters:
        path = cwd
        name = name of the clustered df
    Outputs:
        Saves one csv per ROI to the /dfCreated/updatedCsvs/ folder
    '''

    import pandas as pd

    #get clustered df; contains all cells from ALL ROIs
    dfClust = pd.read_csv(path+'/results/dfCreated/'+name+'.csv', index_col=0)

    #turn off pandas warning for adding values to specified column with 'loc' command
    #per: https://stackoverflow.com/questions/12555323/adding-new-column-to-existing-dataframe-in-python-pandas
    pd.options.mode.chained_assignment = None 

    dfDict = {} #empty dict to store roi:filtered dataframes with new cluster column

    #first break df into different dataframes based on ROI
    #only care about 2 ROIs for this analysis - D16_BB2014A_ROI01 & M27_TT1120A_ROI02
    for roi in ['D16_BB2014A_ROI01','M27_TT1120A_ROI02']:
        #subset to 1 ROI
        dfROI = dfClust[(dfClust['file'] == roi)]

        #get original df based on dfClust['file'] and then iloc the tumor cell using dfClust['index']
        df = pd.read_csv(path+'/data/mIHC_files/'+roi+'.csv', index_col=0)

        #get index and cluster number from dfROI and add these pairs as a tuple to an overall list
        idxClustList = list(zip(dfROI['index'], dfROI['cluster']))

        #add cluster column to filt_df with correct cluster number using idxClustList
        for tup in idxClustList: #for each tuple in idxClustList
            #use loc index and convert to loc's label (name)
            df.loc[tup[0],'cluster'] = tup[1] #tup[0] = original df's .loc index; tup[1]=cluster value of the seed cell

        #add roi:filt_df to dfDict
        dfDict[roi] = df

    #loop through each df again and update the cluster column to map back to the cell class if it is not tumor (and thus is not part of a cluster)    
    for roi,df in dfDict.items():
        df['cluster']=df['cluster'].fillna(df['class']) #if cluster value is a np NaN value, then replace it with its value from the class column

    #save newly updated dfs with their cluster column as new csvs
    for roi,df in dfDict.items():
        df.to_csv(path+'/results/dfCreated/updatedCsvs/'+roi+'_cluster_'+name[16:]+'.csv')
    
    

def clusterCountPerROI(path,name):
    '''
    This function takes in a clustered csv and creates two csvs with the counts of each seed cell per ROI.
    One csv includes counts across all ROIs. One csv includes counts averaged across ROIs.
    The avg csv has both raw count and percentage counts for each cluster.

    Input parameters:
        path = cwd
        name = name of file containing clustered neighborhood data to count cluster abundance from, excluding the '.csv'
    Outputs:
        Saves two csvs saved to 'dfCreated' folder
    '''

    import pandas as pd

    #read clustering csv to analyze (eg. dfNeighClusteredH70allk5; it's a csv that has each seed cell clustered)
    df = pd.read_csv(path+'/results/dfCreated/'+name+'.csv', index_col=0)

    dictCounts = {} #empty dict to store raw counts for each cluster per ROI

    #for each roi in the big df
    for roi in df['file'].unique():

        dfSingleROI = df[df['file'] == roi] #df contains cells from only one region

        #raw frequency counts of cluster
        freq = dfSingleROI['cluster'].value_counts()

        #add to dict to store results
        dictCounts[roi] = freq

    #put counts of clusters into a df
    dfClustCounts = pd.DataFrame(dictCounts).T
    dfClustCounts = dfClustCounts.fillna(0)

    # #save dfClustCounts to csv
    dfClustCounts.to_csv(path+'/results/dfCreated/dfClustCounts'+name[16:]+'_all.csv')

    #create new dataframe with averaged numbers for every region within one patient
    dfClustCounts.index = dfClustCounts.index.str.slice(0,-6)

    #calculate averages on a per patient basis - not on a percent average by regions
    dfClustCountsSum = dfClustCounts.groupby(dfClustCounts.index).sum()
    dfClustCountsSum['Total'] = dfClustCountsSum.sum(axis=1)

    for col in dfClustCountsSum.columns[:-1]: #dont include total
        dfClustCountsSum[str(col)+'_%'] = dfClustCountsSum[col]/dfClustCountsSum['Total'] #divide by the sum of only the first 5 columns

    dfClustCountsAvg = dfClustCountsSum.drop(columns=['Total'])
    
    #save avg df to csv
    dfClustCountsAvg.to_csv(path+'/results/dfCreated/dfClustCounts'+name[16:]+'_avg.csv')


    
def fig3():
    '''
    Single cell analysis of NK cells results in distinct phenotypes related to the proximity to tumor cells and HER2 status.
    This function creates figures 3A-D and supplementary figures S5A-C. 
    Input parameters:
        None  
    Outputs:
        Saves plots to 'figures' folder for figures 3A-D and supplementary figures S5A-C
    '''

    import os
    import plotly.express as px
    import pandas as pd
    import numpy as np
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import fdrcorrection

    print('\n\n***FIGURE 3 - NK cell function versus spatial proximity to neoplastic cells and HER2 status***\n')

    #get current working directory
    path = os.getcwd()
    
    #get csvList
    csvList = getCsvList()
    
    #FIGURE 3A
    file = 'D16_BB2014A_ROI01'
    colorDict = {'Tumor Cells':'rgb(153,153,153)','NK Cells':'rgb(231,41,138)'}
    
    #read specific ROI mIHC csv
    df = pd.read_csv(path+'/data/mIHC_files/'+file+'.csv',index_col=0)
    
    #subset to just tumor and NK
    df = df[df['class'].isin(['Tumor cells','CD56- NKP46+ NK','CD56+ NKP46+ NK','CD56+ NKP46- NK'])]
    df['Cell Type'] = np.where(df['class'] == 'Tumor cells','Tumor Cells','NK Cells')
    
    #plot and save
    fig = px.scatter(df,x='Location_Center_X',y='Location_Center_Y',color='Cell Type',category_orders={'Cell Type':['Tumor Cells','NK Cells']},color_discrete_map=colorDict)
    fig.write_image(path+'/results/figures/figure3A.png')


    #FIGURE 3B
    distThresh = 40 #40 px = 20 µm
    #generate csv storing NK cell function and spatial relationship to neoplastic cells
    nkFunTumSpatial(path=path,csvList=csvList,distThresh=distThresh)
    
    #read csv generated and plot
    colorDict = {'close':'rgb(17,165,121)','far':'rgb(127,60,141)'}
    df = pd.read_csv(path+'/results/dfCreated/dfNKFun_TumorSpatial_all40.csv',index_col=0)
    
    #NOTE: manually change D2_TN1233A patient's PD-1 value to NA - staining was off
    #loc value of 42 and 97 correspond to this patient - .at modifies df - NOTE: these row numbers change if you adjust csvList order
    df.at[42,'PD1'] = np.nan
    df.at[97,'PD1'] = np.nan
    
    #plot close vs far
    fig = px.box(df,y=df.columns[:-4],color='Location',range_y=(-2,102),hover_name='Patient',points='all',color_discrete_map=colorDict,labels={'value':'Percent NK Cells Positive','variable':'Functional Marker'})
    fig.write_image(path+'/results/figures/figure3B.png')
    
    #now run stats on this - Mann-Whitney U + MHT correction
    dfClose = df[df['Location'] == 'close']
    dfFar = df[df['Location'] == 'far']
    
    colList = df.columns[:-4]
    pList = []
    
    #mann-whitney u test
    for col in colList:
        pVal = mannwhitneyu(dfClose[col],dfFar[col]).pvalue
        pList.append(pVal)
    
    #then correct for MHT - Benjamini-Hochberg
    mhtList = fdrcorrection(pList,alpha=0.05)[1]
    print('\nMHT-corrected P-values for Figure 3B:')
    for i in range(len(mhtList)):
        col = colList[i]
        mhtP = mhtList[i]
        print(col+': p =',round(mhtP,3)) 


    #FIGURE 3C,D
    #read in file again
    dfFun = pd.read_csv(path+'/results/dfCreated/dfNKFun_TumorSpatial_all40.csv',index_col=0)
    
    colorDict = {'close':'rgb(17,165,121)','far':'rgb(127,60,141)'}
    
    #subset HER2+/-
    dfFunNeg = dfFun[dfFun['HER2'] == 0]
    dfFunPos = dfFun[dfFun['HER2'] == 1]
    
    dfList = [dfFunNeg,dfFunPos] #first is her2-, then her2+
    titleDict = {0:'HER2-',1:'HER2+'} #for axis naming
    figDict = {0:'Figure 3C',1:'Figure 3D'} #for p-value printing
    
    #loop through her2 status - first fig 3C then fig 3D
    for i in range(len(dfList)):
        
        df = dfList[i]
        figNum = figDict[i][-1]
        
        #NOTE: manually change D2_TN1233A patient's PD-1 value to NA - staining was off
        if i == 1: #this patient is a her2+ patient
            #loc value of 42 and 97 correspond to this patient - .at modifies df - NOTE: the row numbers change if csvList ordering differs
            df.at[42,'PD1'] = np.nan
            df.at[97,'PD1'] = np.nan
        
        fig = px.box(df,y=['KI67','TIM3','PD1'],color='Location',width=500,range_y=[-2,102],points='all',color_discrete_map=colorDict,labels={'value':'Percent NK Cells Positive','variable':titleDict[i]})
        fig.write_image(path+'/results/figures/figure3'+figNum+'.png')
    
        #do stats - for her2+ and then her2- separately
        pList = []
        
        for m in ['KI67','TIM3','PD1']:
            
            #subset by location for each marker
            dfClose = df[df['Location'] == 'close'][m]
            dfFar = df[df['Location'] == 'far'][m]
            #mann-whitney u for the 3 markers
            pVal = mannwhitneyu(dfClose,dfFar).pvalue
            pList.append(pVal)
       
        #correct for MHT - Benjamini-Hochberg
        mhtList = fdrcorrection(pList,alpha=0.05)[1]
        print('\nMHT corrected P-values for '+figDict[i]+':')
        mhtKi67 = mhtList[0]
        mhtTim3 = mhtList[1]
        mhtPd1 = mhtList[2]
        print('KI67: p =',round(mhtKi67,3))
        print('TIM3: p =',round(mhtTim3,3))
        print('PD1: p =',round(mhtPd1,3))

    
    #SUPPLEMENTARY FIGURES S5A-C
    #reread file
    dfFun = pd.read_csv(path+'/results/dfCreated/dfNKFun_TumorSpatial_all40.csv',index_col=0)
    dfFun['Cohort'] = dfFun['Patient'].str[0]
    dfFun['Cohort_HER2'] = dfFun['Cohort']+'_'+dfFun['HER2'].astype(str)
    
    figDict = {'D_0':'A','D_1':'B','M_1':'C'}
    
    for ch in ['D_0','D_1','M_1']:
        #subset to that cohort_her2 group
        dfCH = dfFun[dfFun['Cohort_HER2'] == ch]
        figNum = figDict[ch]
        fig = px.bar(dfCH,x='Patient',y='Total NK Cells',color='Location',barmode='stack',color_discrete_map=colorDict)
        fig.write_image(path+'/results/figures/figureS5'+figNum+'.png')

    print("Figures 3A-D and Supplementary Figures S5A-C saved to 'figures' folder.")
    print('Figure 3 complete.')
      
    

def fig4():
    '''
    Single cell analysis of neoplastic PanCK+ epithelial cells illustrates heterogeneity and high HLA class I expression in close proximity to NK cells.
    This function creates figures 4C-E.
    Input parameters:
        None  
    Outputs:
        Saves plots to 'figures' folder for figures 4C-E
    '''
    
    import os
    import plotly.express as px
    import pandas as pd
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import fdrcorrection

    print('\n\n***FIGURE 4 - Neoplastic cell function versus spatial proximity to NK cells and HER2 status***\n')

    #get current working directory
    path = os.getcwd()
    
    #get csvList
    csvList = getCsvList()
    
    #FIGURE 4C
    distThresh = 40 #40 px = 20 µm
    #generate csv storing NK cell function and spatial relationship to neoplastic cells
    tumorFunNKspatial(path=path,csvList=csvList,distThresh=distThresh)

    #read csv generated and plot
    colorDict = {'close':'rgb(17,165,121)','far':'rgb(127,60,141)'}    
    df = pd.read_csv(path+'/results/dfCreated/dfTumorFun_NKspatial_all40.csv',index_col=0)
    
    #plot close vs far
    fig = px.box(df,y=df.columns[:-4],color='Location',range_y=(-2,102),points='all',color_discrete_map=colorDict,labels={'value':'Percent Tumor Cells Positive','variable':'Functional Marker'})
    fig.write_image(path+'/results/figures/figure4C.png')
    
    #now run stats on this - Mann Whitney U
    dfClose = df[df['Location'] == 'close']
    dfFar = df[df['Location'] == 'far']
    
    #reset with each cohort
    pList = []
    colList = df.columns[:-4]
    
    for col in colList:
        pVal = mannwhitneyu(dfClose[col],dfFar[col]).pvalue
        pList.append(pVal)
    
    #correct for MHT - Benjamini-Hochberg
    mhtList = fdrcorrection(pList,alpha=0.05)[1]
    print('\nMHT corrected P-values for Figure 4C:')
    for i in range(len(mhtList)):
        col = colList[i]
        mhtP = mhtList[i]
        print(col+': p =',round(mhtP,3))        
    
    
    #FIGURE 4D, E
    #read file again    
    dfFun = pd.read_csv(path+'/results/dfCreated/dfTumorFun_NKspatial_all40.csv',index_col=0)
    
    colorDict = {'close':'rgb(17,165,121)','far':'rgb(127,60,141)'}
    
    #subset HER2+/-
    dfFunNeg = dfFun[dfFun['HER2'] == 0]
    dfFunPos = dfFun[dfFun['HER2'] == 1]
    
    dfList = [dfFunNeg,dfFunPos] #first is her2-, then her2+
    titleDict = {0:'HER2-',1:'HER2+'} #for axis naming
    figDict = {0:'Figure 4D',1:'Figure 4E'} #for p-value printing
    
    #loop through her2 status
    for i in range(len(dfList)):
        
        df = dfList[i]
        figNum = figDict[i][-1]
        
        fig = px.box(df,y=['HLA1'],color='Location',points='all',range_y=[-2,102],width=300,height=400,color_discrete_map=colorDict,labels={'value':'Percent Tumor Cells Positive','variable':titleDict[i]})
        fig.write_image(path+'/results/figures/figure4'+figNum+'.png')
    
        #do stats - for her2+ and then her2- separately
        #subset by location for each marker
        dfClose = df[df['Location'] == 'close']['HLA1']
        dfFar = df[df['Location'] == 'far']['HLA1']
    
        pVal = mannwhitneyu(dfClose,dfFar).pvalue
    
        #print her2 status, marker, p value
        print('P-value for '+figDict[i]+':')
        print('HLA1: p =',round(pVal,3))

    print("Figures 4C-E saved to 'figures' folder.")
    print('Figure 4 complete.')        
        


def fig5():
    '''
    Cellular neighborhood clustering of NK cells.
    This function creates figures 5B-F and supplementary figures S7A-C. 
    Input parameters:
        None  
    Outputs:
        Saves plots to 'figures' folder for figures 5B-F and supplementary figures S7A-C
    '''

    import os
    import pandas as pd
    import plotly
    import plotly.express as px
    from seaborn import regplot
    from matplotlib import pyplot as plt
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import fdrcorrection

    print('\n\n***FIGURE 5 - NK spatial cell neighborhoods versus HER2 status***\n')

    #get current working directory
    path = os.getcwd()
    
    #get csvList
    csvList = getCsvList()

    #generate neighborhoods of NK cells
    seedList = ['CD56- NKP46+ NK','CD56+ NKP46- NK','CD56+ NKP46+ NK']
    distThresh = 120 #120px = 60µm
    #make neighborhoods
    makeNeighborhoods(path=path,csvList=csvList,seedList=seedList,distThresh=distThresh)    
    
    #run elbow method to determine optimal number of clusters
    #note that results are not shown in manuscript, so need to manually adjust save parameter
    steps = 16    
    file = 'dfNeighborhoodClusterNK120'
    save = False #toggle to True if saving elbow plot is desired
    elbowMethod(path=path,file=file,steps=steps,save=save)
    #note: from elbowMethod analysis, k=5 is determined
    
    #cluster neighborhoods by cell compositions
    k = 5
    clusterNeighborhoods(path=path,file=file,k=k)

    #FIGURE 5B    
    file = 'dfNeighClusteredNK120k5' #120px = 60µm
    
    #read clustered csv generated above
    df = pd.read_csv(path+'/results/dfCreated/'+file+'.csv', index_col=0)
    
    #filter df to not include the index or file columns
    df = df.drop(['index','file'],axis=1)
    
    #groupby cluster column and take the averages of all of the other columns for each group
    dfCluster = df.groupby(['cluster']).mean()
    
    #reorder columns
    colList = ["countCD8 T cells%", "countCD4 T cells%", "countCD56+ NKP46+ NK%",'countCD56- NKP46+ NK%','countCD56+ NKP46- NK%','countCD11B+ DCs%','countCD11B- DCs%','countCD11B- CD68+ cells%','countMyelomonocytic cells%','countMyeloid other%','countOther CD45+ cells%','countTumor cells%']
    dfCluster = dfCluster[colList]
    
    #rename index clusters to start at 1 rather than 0
    #rename clusters to be order we want- 1 = cd8, 2 = cd4, 3 = other cd45, 4 = tumor/immune, 5 = tumor
    dfCluster = dfCluster.rename(index={0:'1',1:'4',2:'2',3:'5',4:'3'}).sort_values('cluster')
    
    #colors for cell types
    palette = dict(zip(colList,plotly.colors.qualitative.Pastel))
    palette['countTumor cells%'] = 'rgb(239,85,59)'
    
    fig = px.bar(dfCluster,y=dfCluster.columns,barmode='stack',labels={'cluster':'Cluster','value':'Fraction Present'},color_discrete_map=palette)
    fig.update_layout(legend_traceorder="reversed")
    fig.write_image(path+'/results/figures/figure5B.png')
    
    #FIGURE 5C
    #generate new mIHC csvs with cluster annotation for all NK cells
    createCsvsWithClusterCol(path=path,name=file)
    
    #visualize cells in scatterplot for 2 ROIs  
    colList = ['Tumor cells','1','2','3','4','5']
    
    #color map
    palette = {'1':'dodgerblue','2':'magenta','3':'limegreen','4':'darkviolet','5':'coral','Tumor cells':'rgb(153,153,153)','CD8 T cells':'rgb(102,197,204)','CD4 T cells':'rgb(246,207,113)'}
    
    def sorter(column):
        """Sort function"""
        colList = ['Tumor cells','CD8 T cells','CD4 T cells','1','2','3','4','5']
        sortDict = {col: order for order, col in enumerate(colList)}
        return column.map(sortDict)
    
    roiList = ['D16_BB2014A_ROI01','M27_TT1120A_ROI02']
    for roi in roiList:
        file = roi+'_cluster_NK120k5'
    
        #read updated csv with cluster number added
        df = pd.read_csv(path+'/results/dfCreated/updatedCsvs/'+file+'.csv',index_col=0)
    
        df['cluster'] = df['cluster'].replace({'0.0':'1','1.0':'4','2.0':'2','3.0':'5','4.0':'3'})
    
        #only show tumor and NK cell clusters
        df = df[df['cluster'].isin(colList)]
        df = df.sort_values(by='cluster', key=sorter)
    
        #plot scatter reconstructions
        fig = px.scatter(df,x='Location_Center_X',y='Location_Center_Y',color='cluster',color_discrete_map=palette,title='Specimen '+roi[0:3])
        fig.update_traces(marker={'size': 7})
        fig.write_image(path+'/results/figures/figure5C_'+roi[0:3]+'.png')
        
    #SUPPLEMENTARY FIGURE S7A
    name = 'dfNeighClusteredNK120k5'
    #calculate counts of each cluster assignment
    clusterCountPerROI(path=path,name=name)    
    
    #get df from saved csv
    df = pd.read_csv(path+'/results/dfCreated/dfClustCountsNK120k5_all.csv',index_col=0)
    
    #manually rename clusters to match ordering set earlier
    df = df.rename(columns={'0':'1','1':'4','2':'2','3':'5','4':'3'})
    df = df[['1','2','3','4','5']] #reorder columns from 1-5
    sums = df.sum(axis=0)
    total = df.to_numpy().sum()
    dfPerc = pd.DataFrame(index=['1','2','3','4','5'])
    dfPerc['Raw'] = sums.values
    dfPerc['Perc'] = dfPerc['Raw']/total*100 #% out of 100s
    
    fig = px.bar(dfPerc,y='Perc',labels={'index':'Cluster','Perc':'Percent of NK Cell Neighborhood Clusters Present'})
    fig.write_image(path+'/results/figures/figureS7A.png')
    
    #FIGURES 5D, E
    #get df from saved csv
    df = pd.read_csv(path+'/results/dfCreated/dfClustCountsNK120k5_avg.csv',index_col=0)
    
    #generate column list to cluster on based on if there is a % in the column name
    dfPerc = df[df.columns[['%' in col for col in list(df.columns)]]]
    
    #rename columns according to manual order
    colDict = {'0_%':'1','1_%':'4','2_%':'2','3_%':'5','4_%':'3'}
    dfPerc = dfPerc.rename(mapper=colDict, axis=1)
    dfPerc = dfPerc[['1','2','3','4','5']] #order columns from 1-5
    
    #sort descending tumor cluster
    dfPercSort = dfPerc.sort_values('5',ascending=False)
    dfPercSort.index = dfPercSort.index.str[0:-8]
    palette = {'1':'rgb(0,206,209)','2':'rgb(255,0,255)','3':'rgb(50,205,50)','4':'rgb(148,0,211)','5':'rgb(255,127,80)'}
    
    #plot
    fig = px.bar(dfPercSort,y=dfPercSort.columns,barmode='stack',color_discrete_map=palette,labels={'index':'Specimen','value':'Fraction Present','variable':'Cluster'})
    fig.update_layout(legend_traceorder="reversed")
    fig.write_image(path+'/results/figures/figure5D.png')
    
    #look at correlation between sum of t cell clusters vs tumor cluster
    dfPerc['1+2'] = dfPerc['1'] + dfPerc['2']
    
    regplot(x=dfPerc['1+2'],y=dfPerc['5'])
    plt.annotate('Correlation = '+str(round(dfPerc.corr().loc['5','1+2'],3)), (0.6,0.9))
    plt.savefig(path+'/results/figures/figure5E.png',format='png')
    plt.close()

    #SUPPLEMENTARY FIGURE S7B
    #get df from saved csv
    df = pd.read_csv(path+'/results/dfCreated/dfClustCountsNK120k5_all.csv',index_col=0)
    
    #add % column
    for col in df.columns:
        df[str(col)+'_%'] = df[col]/df.iloc[:,0:5].sum(axis=1) #divide by the sum of only the first 5 columns
    
    #generate column list to cluster on based on if there is a % in the column name
    dfPerc = df[df.columns[['%' in col for col in list(df.columns)]]]
    
    colDict = {'0_%':'1','1_%':'4','2_%':'2','3_%':'5','4_%':'3'}
    dfPerc = dfPerc.rename(mapper=colDict, axis=1)
    dfPerc = dfPerc[['1','2','3','4','5']] #order columns from 1-5
    
    #sort by descending tumor cluster
    dfPercSort = dfPerc.sort_values('5',ascending=False)
    
    #color map
    palette = {'1':'darkturquoise','2':'magenta','3':'limegreen','4':'darkviolet','5':'coral'}
    
    #plot
    fig = px.bar(dfPercSort,y=dfPercSort.columns[0:5],barmode='stack',color_discrete_map=palette,labels={'index':'ROI','value':'Fraction Present','variable':'Cluster'})
    fig.update_xaxes(showticklabels=False)
    fig.update_layout(legend_traceorder="reversed")
    fig.write_image(path+'/results/figures/figureS7B.png')

    #FIGURE 5F
    #get df from saved csv
    df = pd.read_csv(path+'/results/dfCreated/dfClustCountsNK120k5_avg.csv',index_col=0)
    
    #generate column list to cluster on based on if there is a % in the column name
    df = df[df.columns[['%' in col for col in list(df.columns)]]]
    
    #get clinical df to add her2 status
    dfClin = pd.read_csv(path+'/data/metadata/clinicalData.csv',index_col=0)
    df['HER2'] = dfClin['HER2']
    
    #reorder columns
    colDict = {'0_%':'1','1_%':'4','2_%':'2','3_%':'5','4_%':'3'}
    df = df.rename(mapper=colDict, axis=1)
    df = df[['1','2','3','4','5','HER2']] #order columns from 1-5
    
    #plot
    fig = px.box(df,y=df.columns[:-1],points='all',color='HER2',labels={'variable':'Cluster','value':'Fraction Present'})
    fig.write_image(path+'/results/figures/figure5F.png')

    #test for significance - HER2
    dfP = df[df['HER2'] == 1]
    dfN = df[df['HER2'] == 0]
    
    pHER2List = [] #to store p values for the HER2 analysis; for MHT correction
    colList = list(df.columns[0:5])
    
    #mann whitney u test
    for col in colList:
        pVal = mannwhitneyu(dfP[col],dfN[col]).pvalue
        pHER2List.append(pVal)
    
    #correct for MHT - Benjamini-Hochberg
    mhtList = fdrcorrection(pHER2List,alpha=0.05)[1]
    
    print('\nMHT corrected P-values for Figure 5F:')
    for i in range(len(mhtList)):
        col = colList[i]
        mhtP = mhtList[i]
        print('Cluster',col+': p =',round(mhtP,3))   
    
    
    #SUPPLEMENTARY FIGURE S7C
    #cohort 1 only: her2+ vs her2-
    #get df from saved csv
    df = pd.read_csv(path+'/results/dfCreated/dfClustCountsNK120k5_avg.csv',index_col=0)
    
    #subset to cohort 1 patients only (duke)
    df = df[['D' in i for i in df.index]]
    
    #generate column list to cluster on based on if there is a % in the column name
    df = df[df.columns[['%' in col for col in list(df.columns)]]]
    
    #get clinical df to add her2 status
    dfClin = pd.read_csv(path+'/data/metadata/clinicalData.csv',index_col=0)
    df['HER2'] = dfClin['HER2']
    
    #reorder columns
    colDict = {'0_%':'1','1_%':'4','2_%':'2','3_%':'5','4_%':'3'}
    df = df.rename(mapper=colDict, axis=1)
    df = df[['1','2','3','4','5','HER2']] #order columns from 1-5
    
    #plot
    fig = px.box(df,y=df.columns[:-1],points='all',color='HER2',labels={'variable':'Cluster','value':'Fraction Present'})
    fig.write_image(path+'/results/figures/figureS7C.png')
    
    #test for significance - HER2
    dfP = df[df['HER2'] == 1]
    dfN = df[df['HER2'] == 0]
    
    pHER2List = [] #to store p values for the HER2 analysis; for MHT correction
    colList = list(df.columns[0:5])
    
    #mann-whitney u
    for col in colList:
        pVal = mannwhitneyu(dfP[col],dfN[col]).pvalue
        pHER2List.append(pVal)
    
    #correct for MHT - Benjamini-Hochberg
    mhtList = fdrcorrection(pHER2List,alpha=0.05)[1]
    
    print('\nMHT corrected P-values for Supplementary Figure S7C:')
    for i in range(len(mhtList)):
        col = colList[i]
        mhtP = mhtList[i]
        print('Cluster',col+': p =',round(mhtP,3))
        
    print("Figures 5B-F and Supplementary Figures S7A-C saved to 'figures' folder.")
    print('Figure 5 complete.')
    print('\nALL ANALYSES COMPLETE.\n')



if __name__=="__main__":
    fig3()
    fig4()
    fig5()
    