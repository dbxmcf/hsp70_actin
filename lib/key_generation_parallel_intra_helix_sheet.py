'''

Code for calculating keys parallely on n-1 CPU cores, where n = Total CPU cores on the system.

@authors: Venkata Sarika Kondra
'''
import os,math,ntpath,socket,argparse
import glob
import time
import pandas as pd
from joblib import Parallel, delayed, cpu_count


parser = argparse.ArgumentParser(description='Parallel Key Generation.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', default='t', \
    help='Name of the sample on which this script should be run.')
parser.add_argument('--thetaBounds', '-theta', metavar='thetaBounds', \
        default='0 , 14.1, 20.31,25.49, 29.99, 34.1, 38, 41.7,45.3, 48.61, 51.91, \
        55.19, 58.29, 61.3, 64.39, 67.3, 70.3, 73.11, 81.69, 84.49, 87.29, 90 ', \
        help='Bin Boundaries for Theta.')
parser.add_argument('--distBounds', '-dist', metavar='distBounds', \
    default='2.81, 7.00, 9.00, 11.00, 14.00, 17.4, 24.16, 30.19, 36.37, 44.78, 175.52', \
    help='Bin Boundaries for maxDist.')
parser.add_argument('--skip', '-skip', metavar='skip', default=False, \
    help='To get only amino acids count make this True')
parser.add_argument('--max_helices', '-max_helices', metavar='max_helices', \
    default=90, help='Helix number in which intra keys are to be calculated')
parser.add_argument('--is_sheet', '-is_sheet', metavar='is_sheet', default=False, \
    help='The type you want to calculate, sheet of helix?')

#Changed according to Sarika's binning
def thetaClass_( binBoundaries, value, type):
        classL = -1
        for i in binBoundaries:
            if value < binBoundaries[0]:#Bins are seperately handled for theta and maxdist.
                if type == 0: #Thetas are not allowed to be less than zero.
                    print(value,binBoundaries[0], 'out of index',binBoundaries )       
                else:
                    classL = binBoundaries.index(binBoundaries[0])+1
                break
            if (value < i) :#If the value is less than the boundary it falls in previous bin.
                if type ==0: classL = binBoundaries.index(i) 
                else: classL =binBoundaries.index(i) + 1
                break
        if value >= binBoundaries[-1]:
            if type ==0:
                if value == binBoundaries[-1]: classL = binBoundaries.index(binBoundaries[-1])
            else : classL = binBoundaries.index(binBoundaries[-1]) +2
        return classL

def calcDist(indexLabel1,indexLabel2):
        x1=xCord[indexLabel1]
        x2=xCord[indexLabel2]
        y1=yCord[indexLabel1]
        y2=yCord[indexLabel2]
        z1=zCord[indexLabel1]
        z2=zCord[indexLabel2]
        distance=(((x1-x2)**2+(y2-y1)**2+(z2-z1)**2)**0.5)
        return distance

def indexFind(index_of_2,i1,j1,k1):
        if index_of_2==i1:
            indexOf0=j1
            indexOf1=k1
        elif index_of_2==j1:
            indexOf0=i1
            indexOf1=k1
        elif index_of_2==k1:
            indexOf0=i1
            indexOf1=j1

        return indexOf0, indexOf1

def processFiles(fileName, max_helices, is_sheet):
    print( fileName)
    allkeys = {}
    all_aminoacidCount = 0
    outFile2 = open((outFolder+"//" +ntpath.basename(fileName)[:-3] + \
            "keys_theta"+str(dTheta)+"_dist"+str(dLen) ), "w")
    # fileTriplets = open((outFolder +"//"+ntpath.basename(fileName)[:-3] + \
    #     "triplets_theta"+str(dTheta)+"_dist"+str(dLen)), "w")
    for helix_in_test in range(1,max_helices):
        #print 'helix: ', helix_in_test
        start_time=time.time()
        filesDict={}
        req_lines=[]
        inFile=open(fileName,'r')
                
        global xCord, yCord, zCord
        aminoAcidName={}
        xCord={}
        yCord={}
        zCord={}
        seq_number={}
        counter=0
        aminoacidCount = []
        helix_count = 1
        helix_dict = {}

        for i in inFile:
            if (i[0:6].rstrip()=="NUMMDL"):
                numOfModels=i[10:14].rstrip()
            if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER')):
                break
            if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
                break
            # For HELIX
            if not is_sheet:
                if df_dict[ntpath.basename(fileName)[:-4].upper()] == 'alpha/alpha barrel':
                    if (i[0:5].rstrip())=="HELIX"  and i[19]=='A':
                        req_lines +=range(int(i.strip()[22:25].strip()),int(i.strip()[34:38].strip())+1)
                        helix_dict[helix_count] = range(int(i.strip()[22:25].strip()),int(i.strip()[34:38].strip())+1)
                        helix_count +=1
            # For SHEET
            if is_sheet:
                if df_dict[ntpath.basename(fileName)[:-4].upper()] == 'beta barrel':
                    if (i[0:5].rstrip())=="SHEET"  and i[32]=='A':
                            req_lines +=range(int(i.strip()[23:26].strip()),int(i.strip()[34:37].strip())+1)
                            helix_dict[helix_count] = range(int(i.strip()[23:26].strip()),int(i.strip()[34:37].strip())+1)
                            helix_count +=1
            
            if(i[0:4].rstrip())=="ATOM" and (i[13:15].rstrip())=="CA" and (i[16]=='A'or i[16]==' ')and i[17:20]!= "UNK" :
                if helix_in_test in helix_dict:
                    if int(i[23:26].rstrip()) in helix_dict[helix_in_test]:
                        #print int(i[23:26].rstrip())
                        aminoacidCount.append(i[17:20])
                        aminoAcidName[counter]=int(aminoAcidLabel[i[17:20]])
                        xCord[counter]=(float(i[30:38]))
                        yCord[counter]=(float(i[38:46]))
                        zCord[counter]=(float(i[46:54]))
                        seq_number[counter]=str(i[22:27])
                        counter+=1
                else:
                    for value_ in allkeys:
                        outFile2.writelines([str(value_),'\t', str(allkeys[value_]),'\n'])
                    print (ntpath.basename(fileName)[:-4],all_aminoacidCount,len(allkeys))
                    return (ntpath.basename(fileName)[:-4],all_aminoacidCount,len(allkeys))
        
        #print 'curl_count:' , helix_dict[helix_in_test]
        print fileName, ':', helix_count
        protLen=len(yCord)
        all_aminoacidCount += protLen
        if not args.skip:
            initialLabel=[]
            sortedLabel=[]
            sortedIndex=[]
            outDist={}
            for m in range(0,3):
                initialLabel.append(0)
                sortedLabel.append(0)
                sortedIndex.append(0)
            for i in range(0,protLen-2):
                for j in range(i+1,protLen-1):
                    for k in range(j+1, protLen):
                        global i1,j1,k1
                        i1=i
                        j1=j
                        k1=k
                        keepLabelIndex={}
                        keepLabelIndex[aminoAcidName[i]]=i
                        keepLabelIndex[aminoAcidName[j]]=j
                        keepLabelIndex[aminoAcidName[k]]=k
                        initialLabel[0]=aminoAcidName[i]
                        initialLabel[1]=aminoAcidName[j]
                        initialLabel[2]=aminoAcidName[k]
                        sortedLabel=list(initialLabel)
                        sortedLabel.sort(reverse=True)
                        if (sortedLabel[0]==sortedLabel[1]) and (sortedLabel[1]==sortedLabel[2]):
                            dist1_2Temp=calcDist(i,j)
                            dist1_3Temp=calcDist(i,k)
                            dist2_3Temp=calcDist(j,k)
                            if dist1_2Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                                indexOf0=i
                                indexOf1=j
                                indexOf2=k
                            elif dist1_3Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                                indexOf0=i
                                indexOf1=k
                                indexOf2=j
                            else:
                                indexOf0=j
                                indexOf1=k
                                indexOf2=i
                        elif(aminoAcidName[i]!=aminoAcidName[j]) and (aminoAcidName[i]!=aminoAcidName[k]) and (aminoAcidName[j]!=aminoAcidName[k]):
                            for index_ in range(0,3):
                                sortedIndex[index_]=keepLabelIndex[sortedLabel[index_]]
                            indexOf0=sortedIndex[0]
                            indexOf1=sortedIndex[1]
                            indexOf2=sortedIndex[2]
                        elif(sortedLabel[0]==sortedLabel[1]) and (sortedLabel[1]!=sortedLabel[2]):
                            indexOf2=keepLabelIndex[sortedLabel[2]]
                            indices=indexFind(indexOf2,i,j,k)
                            a=indexOf2
                            b=indices[0]
                            c=indices[1]
                            dist1_3Temp=calcDist(b,a)
                            dist2_3Temp=calcDist(c,a)
                            if dist1_3Temp>=dist2_3Temp:
                                indexOf0=indices[0]
                                indexOf1=indices[1] 
                            else:
                                indexOf0=indices[1]
                                indexOf1=indices[0]
                        elif(sortedLabel[0]!=sortedLabel[1]) and (sortedLabel[1]==sortedLabel[2]):
                            indexOf0=keepLabelIndex[sortedLabel[0]]
                            indices=indexFind(indexOf0,i,j,k)
                            if calcDist(indexOf0,indices[0])>= calcDist(indexOf0,indices[1]):
                                indexOf1=indices[0]
                                indexOf2=indices[1] 
                            else:
                                indexOf2=indices[0]
                                indexOf1=indices[1]
                        dist01=calcDist(indexOf0,indexOf1)
                        s2=dist01/2
                        dist02=calcDist(indexOf0,indexOf2)
                        s1=dist02
                        dist12=dist01
                        dist03=calcDist(indexOf1,indexOf2)
                        maxDist=max(dist01,dist02,dist03)
                        s3=(((xCord[indexOf0]+xCord[indexOf1])/2-xCord[indexOf2])**2+((yCord[indexOf0]+yCord[indexOf1])/2-yCord[indexOf2])**2+((zCord[indexOf0]+zCord[indexOf1])/2-zCord[indexOf2])**2)**0.5
                        Theta1=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14
                        if Theta1<=90:
                            Theta=Theta1
                        else:
                            Theta=abs(180-Theta1)
                        classT1=thetaClass_(thetaBounds,Theta,0)
                        classL1=thetaClass_(distBounds, maxDist,1)

                        ##getting the positions of AminoAcids in sequence
                        position0 = str(seq_number.values()[indexOf0])
                        position1 = str(seq_number.values()[indexOf1])
                        position2 = str(seq_number.values()[indexOf2])

                        aacd0 = aminoAcidLabel.keys()[aminoAcidLabel.values().index(aminoAcidName[indexOf0])]
                        aacd1 = aminoAcidLabel.keys()[aminoAcidLabel.values().index(aminoAcidName[indexOf1])]
                        aacd2 = aminoAcidLabel.keys()[aminoAcidLabel.values().index(aminoAcidName[indexOf2])]

                        x0 = str(xCord.get(indexOf0))
                        y0 = str(yCord.get(indexOf0))
                        z0 = str(zCord.get(indexOf0))

                        x1 = str(xCord.get(indexOf1))
                        y1 = str(yCord.get(indexOf1))
                        z1 = str(zCord.get(indexOf1))

                        x2 = str(xCord.get(indexOf2))
                        y2 = str(yCord.get(indexOf2))
                        z2 = str(zCord.get(indexOf2))

                        key_2=dLen*dTheta*(numOfLabels**2)*(aminoAcidName[indexOf0]-1)+dLen*dTheta*(numOfLabels)*(aminoAcidName[indexOf1]-1)+dLen*dTheta*(aminoAcidName[indexOf2]-1)+dTheta*(classL1-1)+(classT1-1)
                       
                        if key_2 in allkeys:
                            allkeys[key_2]+=1
                        else:
                            allkeys[key_2]=1

                        # line = (str(key_2)+"\t"+str(aacd0)+"\t"+str(position0)+"\t"+str(aacd1)+"\t"+str(position1)+"\t"+str(aacd2)+"\t"+str(position2)+"\t"+str(classT1)+"\t"+str(Theta)+"\t"+str(classL1)+"\t"+str(maxDist)+"\t"+x0+"\t"+y0+"\t"+z0+"\t"+x1+"\t"+y1+"\t"+z1+"\t"+x2+"\t"+y2+"\t"+z2+"\n")
                        # fileTriplets.writelines(line)        

        end_time=time.time()
        total_time=((end_time)-(start_time))
        #print("NUM OF AMINOACIDS IN HELIX=",protLen)
    for value_ in allkeys:
        outFile2.writelines([str(value_),'\t', str(allkeys[value_]),'\n'])
    print (ntpath.basename(fileName)[:-4],all_aminoacidCount,len(allkeys))
    return (ntpath.basename(fileName)[:-4],all_aminoacidCount,len(allkeys))


if __name__ == '__main__':
    """Executable code starts here."""

    args = parser.parse_args()

    thetaBounds = list(map(float, args.thetaBounds.split(',')))
    dTheta = len(thetaBounds) - 1
    distBounds = list(map(float, args.distBounds.split(',')))
    dLen = len(distBounds) + 1
    numOfLabels = 20
    if socket.gethostname().startswith('qb'):
        path = '/work/sarikak/TSR/Protien_Database/'
        #path = '/work/wxx6941/TSR/'
    else:
        path = '/home/linc/c00219805/Research/Protien_Database/'
    subfolder= '/extracted_new_samples/testing/'+args.sample_name+'//'
    #subfolder= '//'+args.sample_name+'//'
    fileType = "*.pdb"  #"*.ent"#       
    aminoAcidCode=open(path+"aminoAcidCode_lexicographic _new.txt","r") 
    aminoAcidLabel={}
    for amino in aminoAcidCode:
        amino=amino.split()
        aminoAcidLabel[amino[0]]=int(amino[1])
    aminoAcidCode.close()
    outFolder = path+subfolder+"theta"+str(dTheta)+"_dist"+str(dLen)
    #Create output directory
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)   

    df1 = pd.read_csv(path+subfolder+'sample_details.csv')
    df_dict = dict(zip(df1.protein,df1.group))

    print("Working on " + "theta "+str(dTheta)+" dist "+str(dLen))
    files=glob.glob(path+subfolder+fileType)    
    result = Parallel(n_jobs=cpu_count() - 1, verbose=10, backend="multiprocessing", batch_size="auto")(delayed(processFiles)(fileName,args.max_helices, args.is_sheet) for fileName in files)
    
    df2 = pd.DataFrame(result,columns = ['protein','#amino acids','#keys'])
    df_merge = df1.merge(df2,on = 'protein')
    if args.is_sheet:
        df_merge.to_csv(path+subfolder+'sample_details_intra_sheet.csv')
    else:
        df_merge.to_csv(path+subfolder+'sample_details_intra_helix.csv')
    
    #df2.to_csv(path+subfolder+'sample_details_intra_helix.csv')

    print("Parallel Key Generation completed.")
    