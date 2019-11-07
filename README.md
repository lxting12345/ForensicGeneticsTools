# ForensicGeneticsTools
This package which is writed in perl is used to caculate paternity index(PI) in forensic genetics. <br>
it can be used in duo case(father-child, or mother-child), tro case(mother-child-father).

## Usage: <br>
perl PI_v9.pl -table path/to/allelicFrequencyTable -in path/to/genofile -out path/to/resultfile [option].<br>
for example:<br>
perl PI_v9.pl -table ./AF/CNAS_AF_39.txt -in ./test_data/2019-010.txt -out 2019-010-v9.result

required parameter: <br>

|      |      |
|------|------|
|-table|allelic frequency table file| 
|-in|genotype file|
|-out|result file|

options:<br>

|               |            |
|---------------|------------|
|-father_mutation|father_mutation, default is 0.002| 
|-mother_mutation|mother mutation, default is 0.002 in duo, and 0.0005 in tro|
|-precision|accuracy of display, default is 4 digital|
|-float|display in float|
|-scientific|display in sicentific, dufault is scientific|


## format specification <br>
### 1.allelic frequency table <br>
the first row is #name. from the second row, the first colum is name of locus, and allele frequency, each column is seperated by "\t".for example：<br>
#CNAS_AF39 <br>
D3S1358	11=0.0005	12=0.0014	13=0.0014	14=0.0473	15=0.3453	16=0.3277	17=0.2062	18=0.0636	19=0.0061	20=0.0005 <br>
vWA	13=0.0020	14=0.2567	15=0.0303	16=0.1644	17=0.2361	18=0.1947	19=0.0948	20=0.0192	21=0.0018 <br>
FGA	16=0.0002	17=0.0016	18=0.0181	19=0.0445	20=0.0458	21=0.1072	21.2=0.0032	22=0.1866	22.2=0.0038	23=0.2237	23.2=0.0102	<br>

in the AF dir, there is a file named CNAS_AF_39.txt for example.
you can used this table, or also you can build allele frequency table by yourself. <br>

### 2.genotype file <br>
it contain genotype of father and/or mother and child. first col is name of locus, second col is genotype of mother, third col
is child and last col is father, each col is seperated by "\t". if there is no father genotype, use -/- instead. 
it is same to mother. for example: <br>
D3S1358	15/15	15/15	-/-		<br>		
D13S317	8/8	8/8	-/-				<br>
D7S820	10/10	10/10	-/-   <br>

## algorithm description <br>
The principle is according to GB/T37223-2018 <br>
### in duo(二联体)：<br>
PI = (mother pribability*father probability) / (random femal* random male) <br>
### in tro(三联体): <br>
it contain 3 case:
1.you have knowned that the mother(alleged mother) is biological mother, 
  and you want to know wether the father(alleged father) is biological father.<br>
2.you have knowned that the father(alleged father) is biological father, 
  and you want to know wether the mother(alleged mother).<br>
3.you don't know wether the moher and father is biological, 
  and you want to know both. <br>
  
so in tro case, we first judge mother-child as duo, if The PI_MC is greater than 10000, 
we can judge that the mother is biological mother.<br> 
then PI = (motehr probility * father probility)/(mother probility * random male) <br>

otherwise we take father-child as duo, if PI_FC is greater than 10000,
we can judge that the father is biological father. <br> then 
PI = (mother probility * father probility)/(random femal * father probility) <br>

if PI_FC is smaller than 10000. we can't caculate PI_tri, and we only give you the PI of mother-child and PI of father-child.
 
## future
1. finish the perl/tk
2. change code into python3
3. finish the GUI in python3



