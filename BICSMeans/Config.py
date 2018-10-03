#! /usr/bin/env python
import numpy as np
import os, sys
import re, string
from optparse import OptionParser
import tempfile

APP_DIR=os.path.split(os.path.abspath(__file__))[0]
EXE_PATH=os.path.join(APP_DIR, 'source', 'adaptive_bic_0.67')
EXE_PATH2=os.path.join(APP_DIR, 'source', 'SKmeans.0.18')
sys.path.append(os.path.join(APP_DIR,'source'))



class Config():
    def __init__( self, argv ):
        usage='usage: BicSKmeans_0.2 [options]'
        parser = OptionParser(usage="%prog[-f][-q]", version='%prog0.2')
        parser.add_option('-f', dest='file', type='string', help='input file name, Format(rowNames: gene; colNames: array)')
        parser.add_option('--kg', dest='kg', type='int', default='1', help='clusterNum of row. Bic is needed(default, or kg=1); no gene clustering of row(kg=0); the number of clusters of the row for SKmeans(kg>1)')
        parser.add_option('--ka', dest='ka', type='int', default='1', help='clusterNum of column. Bic is needed(default, or ka=1); no sample clustering of column(ka=0); the number of clusters of the row for SKmeans(ka>1)')
        parser.add_option('--cg', dest='cg', type='int', default=50, help='the max cluster number for genes (default 50), which should be a multiple of 5 and less than the number of genes')
        parser.add_option('--ca', dest='ca', type='int', default=50, help='the max cluster number for samples (default 50), which should be a multiple of 5 and less than the number of samples')
        parser.add_option('-r', dest='repeat', type='int', default=100, help='repeat N times for kmeans to find best clusters(default 100)')
        parser.add_option('--eg', dest='eg', type='float', default=1.0, help='reg factor for genes (default 1.0)')
        parser.add_option('--es', dest='es', type='float', default=1.0, help='reg factor for samples (default 1.0)')
        parser.add_option('-l', action='store_true', default=False, help='Specifies to log-transform the data before clustering(default is no log-transform)')
        parser.add_option('--ng', action='store_true', default=False, help='Specifies to normalize each row (gene) in the data (default is no normalization)')
        parser.add_option('--na', action='store_true', default=False, help='Specifies to normalize each column (sample) in the data (default is no normalization)')
        parser.add_option("--hc", action="store_true", default=False, help='Specifies to do hierarchical cluster for clusters from Super Kmeans (default, no hierarchical cluster)')


        (self.options, args) = parser.parse_args(argv)
        # error: Must provide the input file name!
        if self.options.file == None:
            parser.error('-h for help or provide the input file name!')
        if self.options.kg == 0 and self.options.ka == 0:
            parser.error('Please provide kg or ka or both!')

        self.fnConfig=''
        self.data=''
        self.cNames=''
        self.rNames=''
        self.rNames_list=''
        self.data_g=[]
        self.data_s=[]
        self.rNames_s=''
        self.cNames_s=''
        self.cluster_number_g=''
        self.cluster_number_s=''
        self.write()

    def write( self ):

        file=self.options.file
        self.fnConfig_g=os.path.splitext(file)[0]+'_Bic_parameters_g.txt'
        self.fnConfig_s=os.path.splitext(file)[0]+'_Bic_parameters_s.txt'
        self.fnMatrix=os.path.splitext(self.options.file)[0]+'.matrix'
        self.fnMatrix_g=os.path.splitext(self.options.file)[0]+'_g.matrix'
        self.fnMatrix_s=os.path.splitext(self.options.file)[0]+'_s.matrix'
        self.fnClusters_g=os.path.splitext(file)[0]+'_Bic.clusters_g'
        self.fnClusters_s=os.path.splitext(file)[0]+'_Bic.clusters_s'
        self.fnAssignment_g=os.path.splitext(file)[0]+'_Bic.assignment_g'
        self.fnAssignment_s=os.path.splitext(file)[0]+'_Bic.assignment_s'
   
        data=Array(file)
        self.data=data
        self.cNames=data.cNames
        self.rNames=data.rNames
        self.rNames_list=data.header_s

        # define matrix files
        fMat=open(self.fnMatrix,'w')
        fMat_g=open(self.fnMatrix_g, 'w')
        fMat_s=open(self.fnMatrix_s, 'w')
        # l: log
        self.data_g=data.data
        if self.options.l == True:
            if data.forl == 1:
                print 'Data can not be log-transformed!'
                sys.exit()
            else:
                self.data_g=np.log(self.data_g)
        # ka  
#        if self.options.ka != 0:
        self.data_s=self.data_g.transpose()
        self.rNames_s=self.cNames
           # print self.rNames_s
        self.cNames_s=self.rNames
                
        # ng and na
        if self.options.ng == True:
            data_ng=[]
            for l in self.data_g:
                mean = np.mean(l)
                std = np.std(l, ddof=1)
                data_ng.append([float((float(i)-mean)/std) for i in l])
            self.data_g=np.array(list(data_ng))
            self.data_s=self.data_g.transpose()
        if self.options.na == True:
            data_na=[]
            for l in self.data_g.transpose():
                mean = np.mean(l)
                std = np.std(l, ddof=1)
                data_na.append([float((float(i)-mean)/std) for i in l])
            self.data_s=np.array(list(data_na))
############# write fMat file########################################################
        fMat.write('\t'.join([str(i) for i in data.header_g])+'\n')
        line = 0
        #print('data.header_g:'+data.header_g)
        for l in self.data_g:
            fMat.write(data.header_s[line]+'\t'+'\t'.join([str(i) for i in l])+'\n')
            line+=1
        fMat.close()
########################################################################################
############## write fMat_g and fMat_s file ############################################
        for l in self.data_g:
            fMat_g.write('\t'.join([str(i) for i in l])+'\n')
        fMat_g.close()
        if self.data_s != None:
            for l in self.data_s:
                fMat_s.write('\t'.join([str(i) for i in l])+'\n')
            fMat_s.close()
#        else:
#            os.system('rm '+ fMat_s)
  #      if self.options.kg == 0:
   #         os.system('rm '+ fMat_g)
                    
#########################################################################################        
            #write config file
        fg=open(self.fnConfig_g, 'w')
        parameters_g=''
        fs=open(self.fnConfig_s, 'w')
        parameters_s=''

            #file name
        parameters_g+= 'file_name\n'+self.fnMatrix_g+'\n'
        parameters_s+= 'file_name\n'+self.fnMatrix_s+'\n'

            #Input Num
        parameters_g+= 'Input_Num\n'+str(data.nRow)+'\n'
        parameters_s+= 'Input_Num\n'+str(data.nCol)+'\n'

            #Max_Cluster_Num
        if int(self.options.cg) > 85:
            parameters_g+= 'Max_Cluster_Num\n'+str(int(self.options.cg)+15)+'\n'
        else:
            parameters_g+= 'Max_Cluster_Num\n100\n'
        if int(self.options.ca) > 85:
            parameters_s+= 'Max_Cluster_Num\n'+str(int(self.options.ca)+15)+'\n'
        else:
            parameters_s+= 'Max_Cluster_Num\n100\n'

            #Attributes
        parameters_g+= 'Attributes\n'+str(data.nCol)+'\n'
        parameters_s+= 'Attributes\n'+str(data.nRow)+'\n'

            #Max_iteration
        parameters_g+= 'Max_iteration\n40000\n'
        parameters_s+= 'Max_iteration\n40000\n'

            #Trial_num_A
        parameters_g+= 'Trial_num_A\n30\n'
        parameters_s+= 'Trial_num_A\n30\n'
 
            #Trial_num_B
        parameters_g+= 'Trial_num_B\n50\n'
        parameters_s+= 'Trial_num_B\n50\n'

            #is_select
        parameters_g+= 'is_select\n0\n'
        parameters_s+= 'is_select\n0\n'

            #SelectNum
        parameters_g+= 'SelectNum\n1500\n'
        parameters_s+= 'SelectNum\n1500\n'

            #Selected_IDs
        parameters_g+= 'Selected_IDs\ngeneid_g.txt\n'
        parameters_s+= 'Selected_IDs\nsampleid_s.txt\n'

            #is_preprocess
        parameters_g+= 'is_preprocess\n0\n'
        parameters_s+= 'is_preprocess\n0\n'

            #cluster_name
        parameters_g+= 'cluster_name\n'+self.fnClusters_g+'\n'
        parameters_s+= 'cluster_name\n'+self.fnClusters_s+'\n'

            #assignment_name
        parameters_g+= 'assignment_name\n'+self.fnAssignment_g+'\n'
        parameters_s+= 'assignment_name\n'+self.fnAssignment_s+'\n'

            #is_model_selection
        parameters_g+= 'is_model_selection\n1\n'
        parameters_s+= 'is_model_selection\n1\n'
        
            #Cluster_Number
        parameters_g+= 'Cluster_Number\n'+str(self.options.kg)+'\n'
        parameters_s+= 'Cluster_Number\n'+str(self.options.ka)+'\n'
        
            #early_stopping
        parameters_g+= 'early_stopping\n0\n'
        parameters_s+= 'early_stopping\n0\n'

            #min_gap
        parameters_g+= 'min_gap\n2\n'
        parameters_s+= 'min_gap\n2\n'

            #scan_gap_const
        parameters_g+= 'scan_gap_const\n5\n'
        parameters_s+= 'scan_gap_const\n5\n'

            #max_test_num
        if int(self.options.cg) < int(data.nRow):
            parameters_g+= 'max_test_num\n'+str(int((self.options.cg+5)/5))+'\n'
        else:
            parameters_g+= 'max_test_num\n'+str(int(data.nRow/5))+'\n'
        if int(self.options.ca) < int(data.nCol):
            parameters_s+= 'max_test_num\n'+str(int((self.options.ca+5)/5))+'\n'
        else:
            parameters_s+= 'max_test_num\n'+str(int(data.nCol/5))+'\n'

            #reg_factor
        parameters_g+= 'reg_factor\n'+str(self.options.eg)+'\n'
        parameters_s+= 'reg_factor\n'+str(self.options.es)+'\n'

         
        fg.write(parameters_g)
        fg.close()
        fs.write(parameters_s)
        fs.close() 

    def run( self ):
        file = self.options.file
        data=Array(file)
        
#        if self.options.kg != 1:
#            cmd= 'rm '+self.fnConfig_g
#            os.system(cmd)
#        if self.options.ka != 1:
#            cmd= 'rm '+self.fnConfig_s
#            os.system(cmd)
        if self.options.kg == 1:
            cmd=EXE_PATH+' '+self.fnConfig_g+' '+'> '+os.path.splitext(file)[0]+'_Bic_result_g.txt'
            print cmd
            os.system(cmd)
            print 'Bic_reuslt_g.txt finished!'
            for line in open(os.path.splitext(file)[0]+'_Bic_result_g.txt'):
                if re.match('true', line):
                    self.cluster_number_g=line.strip().split(';')[0].split('=')[1].strip()
                    print self.cluster_number_g
#                    cmd='rm '+os.path.splitext(file)[0]+'_Bic_result_g.txt'
#                    os.system(cmd)
#            os.system('rm '+self.fnConfig_g)
#            os.system('rm '+self.fnClusters_g)
#            os.system('rm '+self.fnAssignment_g)
                        
        if self.options.ka == 1:
            cmd=EXE_PATH+' '+self.fnConfig_s+' '+'> '+os.path.splitext(file)[0]+'_Bic_result_s.txt'
            print cmd
            os.system(cmd)
            print 'Bic_reuslt_s.txt finished!'
            for line in open(os.path.splitext(file)[0]+'_Bic_result_s.txt'):
                if re.match('true', line):
                    self.cluster_number_s=line.strip().split(';')[0].split('=')[1].strip()
                    print self.cluster_number_s
#                    cmd='rm '+os.path.splitext(file)[0]+'_Bic_result_s.txt'
#                    os.system(cmd)
#            os.system('rm '+self.fnConfig_s)
#            os.system('rm '+self.fnClusters_s)
#            os.system('rm '+self.fnAssignment_s)
            
        
        self.fnConfigSK_g=os.path.splitext(file)[0]+'_SK_parameters_g.txt'
        self.fnConfigSK_s=os.path.splitext(file)[0]+'_SK_parameters_s.txt'
        self.fnClustersSK_g=os.path.splitext(file)[0]+'_SK.clusters_g'
        self.fnClustersSK_s=os.path.splitext(file)[0]+'_SK.clusters_s'
        self.fnAssignmentSK_g=os.path.splitext(file)[0]+'_SK.assignment_g'
        self.fnAssignmentSK_s=os.path.splitext(file)[0]+'_SK.assignment_s'
        
        if self.options.kg != 0:
#################### write SuperKmeans config file for genes! #################################################
            # write config file
            fgSK=open(self.fnConfigSK_g, 'w')
            configSK_g=''

            # file_name
            configSK_g+= 'file_name\n'+self.fnMatrix_g+'\n'

            # cluster_name
            configSK_g+= 'cluster_name\n'+self.fnClustersSK_g+'\n'

            # assignment_name
            configSK_g+= 'assignment_name\n'+self.fnAssignmentSK_g+'\n'
        
            # Input_Num
            configSK_g+= 'Input_Num\n'+str(data.nRow)+'\n'

            # Cluster_Num
            if self.options.kg != 1:
                configSK_g+= 'Cluster_Num\n'+str(self.options.kg)+'\n'
            else:
                configSK_g+= 'Cluster_Num\n'+str(self.cluster_number_g)+'\n'
       
            # Attributes
            configSK_g+= 'Attributes\n'+str(data.nCol)+'\n'

            # Max_iteration
            configSK_g+= 'Max_iteration\n40000\n'
       
            # Trial_num
            configSK_g+= 'Trial_num\n'+str(self.options.repeat)+'\n'

            fgSK.write(configSK_g)
            fgSK.close()
            
            cmd= EXE_PATH2+' '+self.fnConfigSK_g
            print cmd
            os.system(cmd)

#            os.system('rm '+self.fnConfigSK_g)
#            os.system('rm '+self.fnClustersSK_g)
#            os.system('rm '+self.fnMatrix_g)

        if self.options.ka != 0:
###################### Write SuperKmeans config file for samples! ###################################################
            # write config file
            fsSK=open(self.fnConfigSK_s, 'w')
            configSK_s=''

            # file_name
            configSK_s+= 'file_name\n'+self.fnMatrix_s+'\n'

            # cluster_name
            configSK_s+= 'cluster_name\n'+self.fnClustersSK_s+'\n'

            # assignment_name
            configSK_s+= 'assignment_name\n'+self.fnAssignmentSK_s+'\n'
        
            # Input_Num
            configSK_s+= 'Input_Num\n'+str(data.nCol)+'\n'

            # Cluster_Num
            if self.options.ka != 1:
                configSK_s+= 'Cluster_Num\n'+str(self.options.ka)+'\n'
            else:
                configSK_s+= 'Cluster_Num\n'+str(self.cluster_number_s)+'\n'
       
            # Attributes
            configSK_s+= 'Attributes\n'+str(data.nRow)+'\n'

            # Max_iteration
            configSK_s+= 'Max_iteration\n40000\n'
       
            # Trial_num
            configSK_s+= 'Trial_num\n'+str(self.options.repeat)+'\n'

            fsSK.write(configSK_s)
            fsSK.close()

            cmd= EXE_PATH2+' '+self.fnConfigSK_s
            print cmd
            os.system(cmd)

#            os.system('rm '+self.fnConfigSK_s)
#            os.system('rm '+self.fnClustersSK_s)
#            os.system('rm '+self.fnMatrix_s)

        if self.options.hc == True:
            self.do_hcluster_centers()
#        else:
#            os.system('rm '+self.fnMatrix)
####################################################################################################################

    def do_hcluster_centers(self):

        print 'Transforming data to centers.'

        rNames=[]
        cNames=[]
        data=[]
        
        lines=open(self.fnMatrix).read().rstrip().split('\n')
        cNames=lines[0].split('\t')[1:]

        for line in lines[1:]:
            a=line.split('\t')
            rNames.append(a[0])
            data.append([float(i) for i in a[1:]])

        del lines

        data_raw=np.array(data[:])
        data=np.array(data)

        cluster2gene={}
        gene2cluster={}
        cluster2array={}
        array2cluster={}

        if self.options.kg !=0 :
            i=0
            for line in open(self.fnAssignmentSK_g).read().rstrip().split('\n'):
                gene2cluster[i]=line
                if line in cluster2gene:
                    cluster2gene[line].append(i)
                else:
                    cluster2gene[line]=[i]
                i+=1
            for i in cluster2gene:
                me=np.mean(data[cluster2gene[i]],0)
                for j in cluster2gene[i]:
                    data[j]=me

        if self.options.ka != 0:
            i=0
            for line in open(self.fnAssignmentSK_s).read().rstrip().split('\n'):
                array2cluster[i]=line
                if line in cluster2array:
                    cluster2array[line].append(i)
                else:
                    cluster2array[line]=[i]
                i+=1
            for i in cluster2array:
                me=np.mean(data[:,cluster2array[i]],1)
                for j in cluster2array[i]:
                    data[:,j]=me

        line_out=''
        for i in cNames:
            line_out+='\t'+i
        line_out+='\n'

        for i in range(len(rNames)):
            line_out+=str(rNames[i])+'\t'+string.join([str(j) for j in data[i]],'\t')+'\n'

        open('data2centers.txt','w').write(line_out)

        del data

        print 'Performing hierarchical clustering'
        cmd_line="cluster -f data2centers.txt -m c"
        if self.options.kg != 0:
            cmd_line+=" -g 7"
        if self.options.ka != 0:
            cmd_line+=" -e 7"
        os.system(cmd_line)
        os.system('rm data2centers.txt')

        print 'Summarizing.'
        line_out=''
        line_out2='\tGROUP\n'
        fn2=self.fnMatrix+'.hc'
        ngene=0
        hc_gene_cluster=0
        old_array_cluster=''
        old_gene_cluster=''
        for line in open('data2centers.cdt').read().rstrip().split('\n'):
            a=line.split('\t')
            if a[0]=='GID' or (not a[0]):
                line_out+=line+'\n'
            
            elif a[0]=='AID':
                line_out+=line+'\n'
                arrays=[int(i.split('ARRY')[1].split('X')[0]) for i in a[4:]]
                if self.options.ka != 0:
                    line_out3='\tGROUP\n'
                    hc_array_cluster=0
                    for i in arrays:
                        if old_array_cluster!=array2cluster[i]:
                            hc_array_cluster+=1
                        line_out3+=cNames[i]+'\t'+str(hc_array_cluster)+'\n'
                        old_array_cluster=array2cluster[i]
                    open(fn2+'_K_A'+str(len(cluster2array))+'.kag','w').write(line_out3)

            elif a[0]=='EWEIGHT':
                line_out+=line+'\n'

            else:
                if self.options.kg != 0:
                    rowheaderlen = 4
                else:
                    rowheaderlen = 3
                line_out+='\t'.join(a[:int(rowheaderlen)])
                if self.options.kg != 0:
                    try:
                        gene=int(a[0].split('GENE')[1].split('X')[0])
                    except:
                        gene=ngene
                else:
                    gene = ngene
                ngene+=1
                for i in range(len(a[int(rowheaderlen):])):
                    try:
                        arrayID=arrays[i]
                    except:
                        arrayID=i
                    line_out+='\t'+str(data_raw[gene,arrayID])
                line_out+='\n'
                if self.options.kg != 0:
                    if old_gene_cluster!=gene2cluster[gene]:
                        hc_gene_cluster+=1
                    line_out2+=a[1]+'\t'+str(hc_gene_cluster)+'\n'
                    old_gene_cluster=gene2cluster[gene]
        if self.options.kg != 0:
            open(fn2+'_K_G'+str(len(cluster2gene))+'.kgg','w').write(line_out2)
        fn2_k_cdt=fn2+'_K'
        if self.options.kg != 0:
            fn2_k_cdt+='_G'+str(len(cluster2gene))
        if self.options.ka != 0:
            fn2_k_cdt+='_A'+str(len(cluster2array))
        fn2_k_cdt+='.cdt'
        open(fn2_k_cdt,'w').write(line_out)
        open(fn2+'.cdt','w').write(line_out)

        if self.options.ka != 0:
            os.system('mv data2centers.atr '+fn2+'.atr')            
        if self.options.kg != 0:
            os.system('mv data2centers.gtr '+fn2+'.gtr')
        os.system('rm data2centers.cdt')
#        os.system('rm '+self.fnMatrix)
        return 1
######################################################################################################################
class Array():
    def __init__( self, fn, headerFlag=True):
        self.rNames=[]
        self.fn=fn
        self.header=''
        self.forl = 0
        self.data=[]
        self.read(headerFlag)
        self.data=np.array(list(self.data))
        self.header_g=self.header
        self.header=np.array(self.header)
        self.nRow=len(self.data)
        self.nCol=len(self.data[0])
        self.header_s=self.rNames
        self.rNames=np.reshape(self.rNames,(len(self.rNames),1))
        self.cNames=self.header
   
    def read( self, headerFlag ):
        try:
            os.system('dos2unix '+self.fn)
        except:
            pass

        for line in open(self.fn).read().split('\n'):
            blank=0
            words_new=[]
            if headerFlag:
                headerFlag=False
                self.header=line.split('\t')[0:]
                continue
            words=line.split('\t')
            expression=words[1:]
            if len(expression) > 1 and expression.count(expression[0]) == len(expression):
                continue
            else:
                if len(words)<2:
                    continue
                for w in words[1:]:
                    if w == '':
                        blank += 1
                        words_new.append(float(0))
                    else:
                        words_new.append(float(w))
                        if float(w) < 0:
                            self.forl = 1
                if blank >= int(len(words[1:])/2):
                    continue
                elif blank > 0 and blank < int(len(words[1:])/2):
                    mean = np.array(list(words_new)).mean()
                    i=0
                    for w in words[1:]:
                        i+=1
                        if w == '':
                            words[i]=str(mean)
                self.rNames.append(words[0])
                self.data.append([float(i) for i in words[1:]])
        
if __name__ == "__main__":
    Config(sys.argv).run()


