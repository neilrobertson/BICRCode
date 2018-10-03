#! /usr/bin/env python
import numpy as np
import sys, os
import tempfile
from optparse import OptionParser
from Config import Config

APP_DIR=os.path.split(os.path.abspath(__file__))[0]
EXE_PATH=os.path.join(APP_DIR, 'source', 'adaptive_bic_0.67')
EXE_PATH2=os.path.join(APP_DIR, 'source', 'SKmeans.0.18')
sys.path.append(os.path.join(APP_DIR,'source'))


class BicSKmeans():
    def __init__(self, args):
        self.args=args
        self.config=Config(args)
        self.config.run()
        self.writeKgg()

    def writeKgg(self):
        rNames= np.array(list(self.config.rNames_list))
        cNames= np.array(list(self.config.cNames))

        if self.config.options.kg != 0:
            if self.config.options.kg != 1:
                cluster_numberg =self.config.options.kg
                fnKgg= os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.options.kg)+'.kgg'
            else:
                cluster_numberg=self.config.cluster_number_g
                fnKgg= os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.cluster_number_g)+'.kgg'
            clusterIdx_g=np.loadtxt(self.config.fnAssignmentSK_g)
            data= self.config.data_g
            fgKGG=open(fnKgg, 'w')
            fgKGG.write(cNames[0]+'\tGROUP\n')
            
            data_kgg={}
            line=0
            for l in data:
                data_kgg[rNames[line]]='\t'.join([str(i) for i in l])
                line+=1

            result_kgg={}

            ng=0
            for i in clusterIdx_g:
                if i not in result_kgg:
                    result_kgg[i]=[]
                result_kgg[i].append(rNames[ng])
                ng+=1

            if self.config.options.ka == 0:
                if self.config.options.kg != 1:
                    fngCDT= os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.options.kg)+'.cdt'
                else:
                    fngCDT= os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.cluster_number_g)+'.cdt'
                fgCDT=open(fngCDT, 'w')
                if cNames[0] != None:
                    fgCDT.write(cNames[0]+'\tNAME\tGWEIGHT\t'+'\t'.join(cNames[1:])+'\n')
                else:
                    fgCDT.write('GENE\tNAME\tGWEIGHT\t'+'\t'.join(cNames[1:])+'\n')
                fgCDT.write('EWEIGHT\t\t\t'+'\t'.join(['1']*(self.config.data.nCol))+'\n')

                for i in range(0, int(cluster_numberg)):
                    if i in result_kgg:
                        for j in result_kgg[i]:
                            fgKGG.write(j+'\t'+str(i)+'\n')
                            fgCDT.write(j+'\t'+j+'\t1\t'+data_kgg[j]+'\n')                                
                fgCDT.close()
 
                print 'write file:', fnKgg, 'finished!'
                print 'write file:', fngCDT, 'finished!'
                
            else:
                clusterIdx_s=np.loadtxt(self.config.fnAssignmentSK_s)
                if self.config.options.kg != 1:
                    if self.config.options.ka != 1:
                        fnCDT = os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.options.kg)+'_A'+str(self.config.options.ka)+'.cdt'
                    else:
                        fnCDT = os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.options.kg)+'_A'+str(self.config.cluster_number_s)+'.cdt'
                else:
                    if self.config.options.ka != 1:
                        fnCDT = os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.cluster_number_g)+'_A'+str(self.config.options.ka)+'.cdt'
                    else:
                        fnCDT = os.path.splitext(self.config.options.file)[0]+'_K_G'+str(self.config.cluster_number_g)+'_A'+str(self.config.cluster_number_s)+'.cdt'
                fCDT=open(fnCDT, 'w')

                sample_names=[]
                for i in cNames[1:]:
                    sample_names.append(i)
                result_kag={}
                ns=0
                for i in clusterIdx_s:
                    if i not in result_kag:
                        result_kag[i]=[]
                    result_kag[i].append(sample_names[ns])
                    ns+=1
                
                header_new = []
                order=[]
                if self.config.options.ka != 1:
                    cluster_numbers=self.config.options.ka
                else:
                    cluster_numbers=self.config.cluster_number_s
                for i in range(0, int(cluster_numbers)):
                    if i in result_kag:
                        for j in result_kag[i]:
                            header_new.append(j)
                            order.append(sample_names.index(j))
                fCDT.write('GENE\tNAME\tGWEIGHT\t'+'\t'.join([str(d) for d in header_new])+'\n')
                fCDT.write('EWEIGHT\t\t\t'+'\t'.join(['1']*(self.config.data.nCol))+'\n')

                
                for i in range(0, int(cluster_numberg)):
                    if i in result_kgg:
                        for j in result_kgg[i]:
                            fgKGG.write(j+'\t'+str(i)+'\n')
                            values_new = []
                            values=data_kgg[j].split('\t')
                            for o in order:
                                values_new.append(values[o])
                            fCDT.write(j+'\t'+j+'\t1\t'+'\t'.join([str(d) for d in values_new])+'\n')

                fgKGG.close()
                fCDT.close()
 
                print 'write file:', fnKgg, 'finished!'
                print 'write file:', fnCDT, 'finished!'
#            os.system('rm '+self.config.fnAssignmentSK_g)

        if self.config.options.ka != 0 :
            clusterIdx_s=np.loadtxt(self.config.fnAssignmentSK_s)
            datas= self.config.data_s
            if self.config.options.ka != 1:
                cluster_numbers = self.config.options.ka
                fnKag= os.path.splitext(self.config.options.file)[0]+'_K_A'+str(self.config.options.ka)+'.kag'
                if self.config.options.kg == 0:
                    fnsCDT= os.path.splitext(self.config.options.file)[0]+'_K_A'+str(self.config.options.ka)+'.cdt'
            else:
                cluster_numbers = self.config.cluster_number_s
                fnKag= os.path.splitext(self.config.options.file)[0]+'_K_A'+str(self.config.cluster_number_s)+'.kag'
                if self.config.options.kg == 0:
                    fnsCDT= os.path.splitext(self.config.options.file)[0]+'_K_A'+str(self.config.cluster_number_s)+'.cdt'
            if self.config.options.kg == 0:
                fsCDT=open(fnsCDT, 'w')
                fsCDT.write('SAMPLE\tNAME\tGWEIGHT\t'+'\t'.join([str(g) for g in rNames])+'\n')
                fsCDT.write('EWEIGHT\t\t\t'+'\t'.join(['1']*(self.config.data.nRow))+'\n')
            
            fsKAG=open(fnKag, 'w')
            fsKAG.write('SAMPLE\tGROUP\n')

            data_kag={}
            line=0
            for l in datas:
                data_kag[cNames[line+1]]='\t'.join([str(i) for i in l])
                line+=1

            result_kag={}
            sample_names=[]
            for i in cNames[1:]:
                sample_names.append(i)
            ns=0
            for i in clusterIdx_s:
                if i not in result_kag:
                    result_kag[i]=[]
                result_kag[i].append(sample_names[ns])
                ns+=1
            for i in range(0, int(cluster_numbers)):
                if i in result_kag:
                    for j in result_kag[i]:
                        fsKAG.write(j+'\t'+str(i)+'\n')
                        if self.config.options.kg == 0:
                            fsCDT.write(j+'\t'+j+'\t1\t'+data_kag[j]+'\n')                                
            fsKAG.close()
 
            print 'write file:', fnKag, 'finished!'
            if self.config.options.kg == 0:
                fsCDT.close()
                print 'write file:', fnsCDT, 'finished!'
#            os.system('rm '+self.config.fnAssignmentSK_s)
        os.system('rm *.matrix')
        os.system('rm *_g*')
        os.system('rm *_s*')
        try:
            os.system('rm temp_file.txt')
        except:
            pass

        
if __name__ == "__main__":
    BicSKmeans(sys.argv)
