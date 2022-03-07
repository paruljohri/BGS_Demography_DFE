#This is to generate the sampling scheme for f1, f2, and f3:                    
                                                                                
from random import shuffle                                                      
import sys                                                                      
                                                                                
result = open("/home/pjohri1/eqm_disc_parameters+5.txt", 'w+')                  
result.write("simID" + '\t' + "f0" + '\t' + "f1" + '\t' + "f2" + '\t' + "f3" + '\n')
l_paras = []                                                                    
i = 0                                                                           
while i <= 100:                                                                 
    j = 0                                                                       
    while j <= 100:                                                             
        k = 0                                                                   
        while k <= 100:                                                         
            if (i+j+k) <= 100:                                                  
                l_paras.append(str(100-i-j-k) + '\t' + str(i) + '\t' + str(j) + '\t' + str(k) + '\n')
            k = k + 5                                                           
        j = j + 5                                                               
    i = i + 5                                                                   
                                                                                
shuffle(l_paras)
ID = 1                                                                          
for x in l_paras:                                                               
    result.write("sim" + str(ID) + '\t' + x + '\n')                             
    ID = ID + 1                                                                 
result.close()                                                                  
                                                                                
print ("Finished")
