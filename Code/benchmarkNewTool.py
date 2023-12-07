import scipy.stats
import pandas
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import cmocean
import sklearn
import sklearn.metrics


# Parse input parameters
parser = argparse.ArgumentParser(description="Benchmark a new tool against the others in the current database",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--newToolScores",help="Path to file containing variant calls for the tool to evaluate. Should contain columns for 'Chromosome','Position','Ref','Alt', and the tool prediction values")
parser.add_argument("--toolColumnName",help="Header value for the fifth column of the newToolScores file")
parser.add_argument("--toolName",help="Name to be used for the new tool in all outputs")
parser.add_argument("--dataDirectory",help="Path to directory containing the test sets, test sample data, and initial database values")
args = vars(parser.parse_args())

# Establish variables
inputFile=args['newToolScores']
thisToolName=args['toolName']
thisToolColumnName=args['toolColumnName']
workingDir=args['dataDirectory']
standardToolsNameConverter={'Polyphen2_HDIV_score':'Polyphen2 HDIV','Polyphen2_HVAR_score':'Polyphen2 HVAR','MutationTaster_score':'MutationTaster','VEST4_score':'VEST4',
                            'MetaSVM_score':'MetaSVM','MetaLR_score':'MetaLR','M-CAP_score':'M-CAP','REVEL_score':'REVEL','MutPred_score':'MutPred','MVP_score':'MVP','MPC_score':'MPC',
                            'PrimateAI_score':'PrimateAI','DEOGEN2_score':'DEOGEN2','CADD_raw':'CADD','DANN_score':'DANN','fathmm-MKL_coding_score':'fathmm-MKL','fathmm-XF_coding_score':'fathmm-XF',
                            'GenoCanyon_score':'GenoCanyon','integrated_fitCons_score':'fitCons','GERP++_RS':'GERP++','phyloP100way_vertebrate':'phyloP 100-way (vertebrate)',
                            'phyloP470way_mammalian':'phyloP 470-way (mammalian)','phyloP17way_primate':'phyloP 17-way (primate)','phastCons100way_vertebrate':'phastCons 100-way (vertebrate)',
                            'phastCons470way_mammalian':'phastCons 470-way (mammalian)','phastCons17way_primate':'phastCons 17-way (primate)','SiPhy_29way_logOdds':'SiPhy 29-way',
                            'SIFT_score':'SIFT','SIFT4G_score':'SIFT4G','LRT_score':'LRT','FATHMM_score':'FATHMM','PROVEAN_score':'PROVEAN','bStatistic':'bStatistic',
                            'MutationAssessor_score':'MutationAssessor','Eigen-raw_coding':'Eigen','Eigen-PC-raw_coding':'Eigen-PC','Maverick_noZygosity':'MAVERICK (No Zygosity)',
                            'MAPPIN_noZygosity':'MAPPIN (No Zygosity)','BayesDel_addAF_score':'BayesDel with AF','BayesDel_noAF_score':'BayesDel without AF','ClinPred_score':'ClinPred',
                            'LIST-S2_score':'LIST-S2','VARITY_R_score':'VARITY R','VARITY_ER_score':'VARITY ER','VARITY_R_LOO_score':'VARITY R LOO','VARITY_ER_LOO_score':'VARITY ER LOO',
                            'AlphaMissense':'AlphaMissense','gMVP_score':'gMVP','MetaRNN_score':'MetaRNN','ESM1b':'ESM1b','Maverick':'Maverick','MAPPIN':'MAPPIN',
                            'EVE':'EVE','PrimateAI3D':'PrimateAI-3D',thisToolName:thisToolName}

# Create plotting functions   
def plotRanksNorm(normalizedRanks,description='New tool',fileName='tmp.png',maxSize=0.003):
    totalSites=len(normalizedRanks)
    plt.figure(figsize=(15, 10))
    # drop Maverick and Mappin and their NoZygosity versions from the standard set
    standardRanksData=normalizedRanks.drop(columns=['Maverick_noZygosity','MAPPIN_noZygosity','Maverick','MAPPIN'])
    cm=cmocean.cm.phase
    x=np.arange(0,maxSize+0.0001,0.0001)
    yMaverick=np.zeros((len(x)))
    yMappin=np.zeros((len(x)))
    yMaverickNZ=np.zeros((len(x)))
    yMappinNZ=np.zeros((len(x)))
    standardToolsList=standardRanksData.columns.values
    yStandardTools=np.zeros((len(x),len(standardToolsList)))
    for i in range(len(x)):
        yMaverick[i]=100*len(normalizedRanks.loc[((normalizedRanks['Maverick']>0) & (normalizedRanks['Maverick']<=x[i])),:])/len(normalizedRanks.loc[(normalizedRanks['Maverick']>0),:])
        yMappin[i]=100*len(normalizedRanks.loc[((normalizedRanks['MAPPIN']>0) & (normalizedRanks['MAPPIN']<=x[i])),:])/len(normalizedRanks.loc[(normalizedRanks['MAPPIN']>0),:])
        for j in range(len(standardToolsList)):
            yStandardTools[i,j]=100*len(standardRanksData.loc[((standardRanksData[standardToolsList[j]]>0) & (standardRanksData[standardToolsList[j]]<=x[i])),:])/len(standardRanksData.loc[(standardRanksData[standardToolsList[j]]>0),:])
        yMaverickNZ[i]=100*len(normalizedRanks.loc[((normalizedRanks['Maverick_noZygosity']>0) & (normalizedRanks['Maverick_noZygosity']<=x[i])),:])/len(normalizedRanks.loc[(normalizedRanks['Maverick_noZygosity']>0),:])
        yMappinNZ[i]=100*len(normalizedRanks.loc[((normalizedRanks['MAPPIN_noZygosity']>0) & (normalizedRanks['MAPPIN_noZygosity']<=x[i])),:])/len(normalizedRanks.loc[(normalizedRanks['MAPPIN_noZygosity']>0),:])
    plt.plot(x,yMaverick,color="#5960a0", linewidth=5,label='MAVERICK')
    plt.plot(x,yMappin,color="#bd5c3f", linewidth=5,label='MAPPIN')
    plt.plot(x,yMaverickNZ,color="#5960a0",linestyle='dashed', linewidth=5,label='MAVERICK (No Zygosity)')
    plt.plot(x,yMappinNZ,color="#bd5c3f",linestyle='dashed', linewidth=5,label='MAPPIN (No Zygosity)')
    for i in range(len(standardToolsList)):
        plt.plot(x,yStandardTools[:,i],linewidth=2,label=standardToolsNameConverter[standardToolsList[i]],color=cm(float(i)/(len(standardToolsList)-1)))
    plt.ylim([-5,105])
    plt.yticks([0,20,40,60,80,100],fontsize=15)
    plt.title(description,fontsize=25)
    plt.ylabel('Cumulative\npercent of\ncausal\nvariants\nfound',fontsize=20,rotation=0,labelpad=40)
    plt.xlabel('Proportion of variants inspected',fontsize=20)
    plt.legend(bbox_to_anchor=(1.5, 0.5),loc="right",fontsize=10)
    plt.tight_layout()
    plt.savefig(fileName,bbox_inches='tight',dpi=600)
    
def plotAUCScores(normalizedRanks,description='New tool',fileName='tmp.png',maxSize=0.003):
    x=np.arange(0,maxSize+0.0001,0.0001)
    totalSites=len(normalizedRanks)
    standardToolsList=normalizedRanks.columns.values
    y=np.zeros((len(x),len(normalizedRanks)))
    yLabels=[]
    for i in range(len(standardToolsList)):
        yLabels.append(standardToolsNameConverter[standardToolsList[i]])
    for i in range(len(x)):
        for j in range(len(standardToolsList)):
            y[i,j]=len(normalizedRanks.loc[((normalizedRanks[standardToolsList[j]]>0) & (normalizedRanks[standardToolsList[j]]<=x[i])),:])/len(normalizedRanks.loc[(normalizedRanks[standardToolsList[j]]>0),:])
    aucs=np.zeros(len(standardToolsList))
    for i in range(len(standardToolsList)):
        aucs[i]=np.round(sklearn.metrics.auc(x/maxSize,y[:,i]),4)
    
    sortOrder=scipy.stats.rankdata(-aucs,method='ordinal')-1
    aucs2=np.zeros(len(aucs))
    yLabels2=[]
    for i in range(len(aucs)):
        aucs2[i]=aucs[np.where(sortOrder==i)[0][0]]
        yLabels2.append(yLabels[np.where(sortOrder==i)[0][0]])
    plt.figure(figsize=(20, 10))
    x=np.arange(len(aucs2))
    plt.bar(x,aucs2,align='center',color="#5960a0")
    plt.xticks(x,yLabels2,fontsize=10,rotation=90)
    plt.ylim([0,1.05])
    plt.xlim([-0.7,39.7])
    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],fontsize=10)
    plt.ylabel('Area Under\nthe Curve',fontsize=15,rotation=0,labelpad=50)
    plt.title(description,fontsize=25)
    plt.tight_layout()
    plt.savefig(fileName,bbox_inches='tight',dpi=600)

def plotRanksNonNorm(inputRanks,description='Known Genes',fileName='tmp.png',maxSize=100):
    totalSites=len(inputRanks)
    plt.figure(figsize=(15, 10))
    #hsv=plt.get_cmap('hsv')
    # drop Maverick and Mappin and their NZ versions from the standard set
    standardRanksData=inputRanks.drop(columns=['Maverick_noZygosity','MAPPIN_noZygosity','Maverick','MAPPIN'])
    cm=cmocean.cm.phase
    x=np.arange(0,maxSize+1,maxSize/20)
    yMaverick=np.zeros((len(x)))
    yMappin=np.zeros((len(x)))
    yMaverickNZ=np.zeros((len(x)))
    yMappinNZ=np.zeros((len(x)))
    standardToolsList=standardRanksData.columns.values
    yStandardTools=np.zeros((len(x),len(standardToolsList)))
    for i in range(len(x)):
        yMaverick[i]=100*len(inputRanks.loc[((inputRanks['Maverick']>0) & (inputRanks['Maverick']<=x[i])),:])/len(inputRanks.loc[(inputRanks['Maverick']>0),:])
        yMappin[i]=100*len(inputRanks.loc[((inputRanks['MAPPIN']>0) & (inputRanks['MAPPIN']<=x[i])),:])/len(inputRanks.loc[(inputRanks['MAPPIN']>0),:])
        for j in range(len(standardToolsList)):
            yStandardTools[i,j]=100*len(standardRanksData.loc[((standardRanksData[standardToolsList[j]]>0) & (standardRanksData[standardToolsList[j]]<=x[i])),:])/len(standardRanksData.loc[(standardRanksData[standardToolsList[j]]>0),:])
        yMaverickNZ[i]=100*len(inputRanks.loc[((inputRanks['Maverick_noZygosity']>0) & (inputRanks['Maverick_noZygosity']<=x[i])),:])/len(inputRanks.loc[(inputRanks['Maverick_noZygosity']>0),:])
        yMappinNZ[i]=100*len(inputRanks.loc[((inputRanks['MAPPIN_noZygosity']>0) & (inputRanks['MAPPIN_noZygosity']<=x[i])),:])/len(inputRanks.loc[(inputRanks['MAPPIN_noZygosity']>0),:])
    plt.plot(x,yMaverick,color="#5960a0", linewidth=5,label='MAVERICK')
    plt.plot(x,yMappin,color="#bd5c3f", linewidth=5,label='MAPPIN')
    plt.plot(x,yMaverickNZ,color="#5960a0",linestyle='dashed', linewidth=5,label='MAVERICK (No Zygosity)')
    plt.plot(x,yMappinNZ,color="#bd5c3f",linestyle='dashed', linewidth=5,label='MAPPIN (No Zygosity)')
    for i in range(len(standardToolsList)):
        plt.plot(x,yStandardTools[:,i],linewidth=2,label=standardToolsNameConverter[standardToolsList[i]],color=cm(float(i)/(len(standardToolsList)-1)))
    plt.ylim([-5,105])
    plt.yticks([0,20,40,60,80,100],fontsize=15)
    plt.title(description,fontsize=25)
    plt.ylabel('Cumulative\npercent of\ncausal\nvariants\nfound',fontsize=20,rotation=0,labelpad=40)
    plt.xlabel('Proportion of variants inspected',fontsize=20)
    plt.legend(bbox_to_anchor=(1.5, 0.5),loc="right",fontsize=10)
    plt.tight_layout()
    plt.savefig(fileName,bbox_inches='tight',dpi=600)
    
def plotAUCScoresNonNorm(inputRanks,description='Known Genes',fileName='tmp.png',maxSize=20):
    x=np.arange(0,maxSize+1)
    totalSites=len(inputRanks)
    standardToolsList=inputRanks.columns.values
    y=np.zeros((len(x),len(inputRanks)))
    yLabels=[]
    for i in range(len(standardToolsList)):
        yLabels.append(standardToolsNameConverter[standardToolsList[i]])
    for i in range(len(x)):
        for j in range(len(standardToolsList)):
            y[i,j]=len(inputRanks.loc[((inputRanks[standardToolsList[j]]>0) & (inputRanks[standardToolsList[j]]<=x[i])),:])/len(inputRanks.loc[(inputRanks[standardToolsList[j]]>0),:])
    aucs=np.zeros(len(standardToolsList))
    for i in range(len(standardToolsList)):
        aucs[i]=np.round(sklearn.metrics.auc(x/maxSize,y[:,i]),4)
    
    sortOrder=scipy.stats.rankdata(-aucs,method='ordinal')-1
    aucs2=np.zeros(len(aucs))
    yLabels2=[]
    for i in range(len(aucs)):
        aucs2[i]=aucs[np.where(sortOrder==i)[0][0]]
        yLabels2.append(yLabels[np.where(sortOrder==i)[0][0]])
    plt.figure(figsize=(20, 10))
    x=np.arange(len(aucs2))
    plt.bar(x,aucs2,align='center',color="#5960a0")
    plt.xticks(x,yLabels2,fontsize=10,rotation=90)
    plt.ylim([0,1.05])
    plt.xlim([-0.7,39.7])
    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],fontsize=10)
    plt.ylabel('Area Under\nthe Curve',fontsize=15,rotation=0,labelpad=50)
    plt.title(description,fontsize=25)
    plt.tight_layout()
    plt.savefig(fileName,bbox_inches='tight',dpi=600)




#### Main work
variantScores=pandas.read_csv(inputFile,sep='\t',low_memory=False)
# drop non-autosomal variants
variantScores=variantScores.rename(columns={'Chromosome':'hg19_chr','Position':'hg19_pos(1-based)','Ref':'ref','Alt':'alt'})
variantScores=variantScores.loc[((variantScores['hg19_chr']!='chrX') & (variantScores['hg19_chr']!='chrY') & (variantScores['hg19_chr']!='chrMT') & (variantScores['hg19_chr']!='chrM')),:].reset_index(drop=True)
variantScores['hg19_chr']=variantScores.loc[:,'hg19_chr'].str[3:].astype(int)

testSampleVariants=pandas.read_csv(os.path.join(workingDir,'allSamplesExonicVariants.txt'),sep='\t',low_memory=False)
# drop variants that are homozygous ref
testSampleVariants=testSampleVariants.loc[~(testSampleVariants['genotype']==0),:].reset_index(drop=True)
# drop non-autosomal variants
testSampleVariants=testSampleVariants.loc[~(testSampleVariants['hg19_chr'].isin(['X','Y'])),:]
testSampleVariants['hg19_chr']=testSampleVariants.loc[:,'hg19_chr'].astype(int)
sampleIDs=testSampleVariants.loc[:,'sampleID'].drop_duplicates(keep='first').reset_index(drop=True)

variantScores2=variantScores.loc[:,['hg19_chr','hg19_pos(1-based)','ref','alt',thisToolColumnName]]

years=[2017,2018,2019,2020,2021,2022,2023]

for year in years:
    print("Starting year " + str(year))
    rec=pandas.read_csv(os.path.join(workingDir,str(year) + '_autosomalRecessive_notInClinvarAtAllBefore.txt'),sep='\t',low_memory=False)
    rec=rec.rename(columns={'Chromosome':'hg19_chr','Position':'hg19_pos(1-based)','ReferenceAllele':'ref','AlternateAllele':'alt'})
    rec=rec.loc[~(rec['hg19_chr'].isin(['X','Y'])),:]
    rec['hg19_chr']=rec.loc[:,'hg19_chr'].astype(int)
    rec=rec.merge(variantScores2,how='inner',on=['hg19_chr','hg19_pos(1-based)','ref','alt'])
    dom=pandas.read_csv(os.path.join(workingDir,str(year) + '_autosomalDominant_notInClinvarAtAllBefore.txt'),sep='\t',low_memory=False)
    dom=dom.rename(columns={'Chromosome':'hg19_chr','Position':'hg19_pos(1-based)','ReferenceAllele':'ref','AlternateAllele':'alt'})
    dom=dom.loc[~(dom['hg19_chr'].isin(['X','Y'])),:]
    dom['hg19_chr']=dom.loc[:,'hg19_chr'].astype(int)
    dom=dom.merge(variantScores2,how='inner',on=['hg19_chr','hg19_pos(1-based)','ref','alt'])

    # create the score collection variables
    thisToolRanks=pandas.DataFrame(np.zeros((len(sampleIDs)*(len(rec)+len(dom)),2)),columns=[thisToolName,'classLabel'])
    thisToolRanks['classLabel']=1
    thisToolRanksNorm=thisToolRanks.copy(deep=True)

    for i in range(len(sampleIDs)):
        print("Starting Sample #" + str(i) + " of " + str(len(sampleIDs)))
        thisSample=str(sampleIDs[i])
        thisSampleVariants=testSampleVariants.loc[testSampleVariants['sampleID']==thisSample,:].reset_index(drop=True)
        # get the scores for the variants in this sample
        thisSampleVariantsScored=variantScores2.merge(thisSampleVariants,how='inner',on=['hg19_chr','hg19_pos(1-based)','ref','alt'])

        # spike in each recessive gene variant
        thisSampleThisToolVariantsScored=thisSampleVariantsScored.loc[:,thisToolColumnName].astype(float).values
        for j in range(len(rec)):
            thisToolRanks.loc[(i*(len(rec)+len(dom)))+j,thisToolName]=scipy.stats.rankdata(-np.concatenate((float(rec.loc[j,thisToolColumnName]),thisSampleThisToolVariantsScored),axis=None),method='average')[0]
            thisToolRanks.loc[(i*(len(rec)+len(dom)))+j,'classLabel']=2
            thisToolRanksNorm.loc[(i*(len(rec)+len(dom)))+j,thisToolName]=thisToolRanks.loc[(i*(len(rec)+len(dom)))+j,thisToolName]/len(thisSampleVariantsScored)
            thisToolRanksNorm.loc[(i*(len(rec)+len(dom)))+j,'classLabel']=2

        # spike in each dominant gene variant
        thisSampleThisToolVariantsScored=thisSampleVariantsScored.loc[:,thisToolColumnName].astype(float).values
        for j in range(len(dom)):
            thisToolRanks.loc[(i*(len(rec)+len(dom)))+len(rec)+j,thisToolName]=scipy.stats.rankdata(-np.concatenate((float(dom.loc[j,thisToolColumnName]),thisSampleThisToolVariantsScored),axis=None),method='average')[0]
            thisToolRanksNorm.loc[(i*(len(rec)+len(dom)))+len(rec)+j,thisToolName]=thisToolRanks.loc[(i*(len(rec)+len(dom)))+len(rec)+j,thisToolName]/len(thisSampleVariantsScored)

    thisToolRanks.to_csv(os.path.join(workingDir,thisToolName + 'Ranks' + str(year) + 'SpikeIns.txt'),sep='\t',index=False)
    thisToolRanksNorm.to_csv(os.path.join(workingDir,thisToolName + 'Ranks' + str(year) + 'SpikeIns_normalized.txt'),sep='\t',index=False)
    databaseNorm=pandas.read_csv(os.path.join(workingDir,str(year) + '_combinedNormalizedResults.txt.gz'),sep='\t',low_memory=False,compression='gzip')
    databaseNorm[thisToolName]=thisToolRanksNorm.loc[:,thisToolName]
    databaseNorm.to_csv(os.path.join(workingDir,str(year) + '_combinedNormalizedResults.txt.gz'),sep='\t',index=False,compression='gzip')
    databaseNonNorm=pandas.read_csv(os.path.join(workingDir,str(year) + '_combinedResults.txt.gz'),sep='\t',low_memory=False,compression='gzip')
    databaseNonNorm[thisToolName]=thisToolRanks.loc[:,thisToolName]
    databaseNonNorm.to_csv(os.path.join(workingDir,str(year) + '_combinedResults.txt.gz'),sep='\t',index=False,compression='gzip')
    plotRanksNorm(databaseNorm,description="Normalized rank of\ntrue pathogenic variants " + str(year),
                  fileName=os.path.join(workingDir,str(year) + '_normalizedRanks_with' + thisToolName + '.png'),maxSize=0.003)
    plotAUCScores(databaseNorm,description="VarPB Area Under the Curve for " + str(year),
                  fileName=os.path.join(workingDir,str(year) +'_AUC_barplots_normalizedRanks_with' + thisToolName + '.png'),maxSize=0.003)
    plotRanksNonNorm(databaseNonNorm,description="Rank of true pathogenic variants " + str(year),
                     fileName=os.path.join(workingDir,str(year) + '_ranks100_with' + thisToolName + '.png'),maxSize=100)
    plotAUCScoresNonNorm(databaseNonNorm,description="VarPB Area Under the Curve for " + str(year) + " (non-normalized, top 100)",
                         fileName=os.path.join(workingDir,str(year) +'_AUC_barplots_ranks100_with' + thisToolName + '.png'),maxSize=100)
    plotRanksNonNorm(databaseNonNorm,description="Rank of true pathogenic variants " + str(year),
                     fileName=os.path.join(workingDir,str(year) + '_ranks20_with' + thisToolName + '.png'),maxSize=20)
    plotAUCScoresNonNorm(databaseNonNorm,description="VarPB Area Under the Curve for " + str(year) + " (non-normalized, top 20)",
                         fileName=os.path.join(workingDir,str(year) +'_AUC_barplots_ranks20_with' + thisToolName + '.png'),maxSize=20)


