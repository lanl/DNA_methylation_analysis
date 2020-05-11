import numpy as np
import pandas as pd
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.KGML import KGML_pathway
pwy=KGML_parser.read(open("/Users/shounak/Downloads/ko00061.xml",'r'))
for i in demo_pwy.orthologs:
    print (print (demo.name.lstrip("path:")+","+",".join(i.name.split(" ")).replace("ko:",""))
