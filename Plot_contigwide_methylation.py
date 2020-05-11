#!/Users/shounak/anaconda3/bin/python3
import matplotlib.pylab as plb
import matplotlib.colors
import argparse
def write_color_coded_gff(meth_value,cmap):
	return (str(matplotlib.colors.to_hex(cmap(int(255.0*meth_value)))))
gff_cmap=plb.get_cmap("magma") #load the magma color scale we are using in our heatmaps
list_ct=0
optparse=argparse.ArgumentParser()
optparse.add_argument("-c","--csvfile",help="list of methylation ratios")
optparse.add_argument("-g","--gff",help="annotation file to augment")
optparse.add_argument("-o","--outfile",help="output annotation file")
argstr=optparse.parse_args()
gff_file=open(argstr.gff,'r')
gff_out=open(argstr.outfile,'w')
for line in gff_file.readlines():
    gff_out.write(line)
gff_file.close()
meth_list_file=open(argstr.csvfile)
for item in meth_list_file.readlines()[1:]:
    fields=item.rstrip("\n").split(',')
    list_ct+=1
    line_out="\t".join([fields[0],"CRS_SB","methylation_C",fields[1],fields[2],fields[3],fields[4],".","ID="+fields[5]+str(list_ct)+";color="+write_color_coded_gff(float(fields[3]),gff_cmap)+"\n"])
    gff_out.write(line_out)
gff_out.close()