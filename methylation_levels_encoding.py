import pysam,argparse
from multiprocessing import Pool
from subprocess import call

parser = argparse.ArgumentParser(description = "Extracts read-level methylation ratios from Bismark CpG reports")
parser.add_argument("-i","--input",metavar="",type=str,help = "the input bam file which is sorted and indexed ")
parser.add_argument("-r","--region",metavar="",type=str,help =" BED file of genomic regions (coordinates should match reference genome)")
parser.add_argument("-b","--ob",metavar="",type=str,help = "CpG_OB* file (bottom strand) from bismark methylation extractor")
parser.add_argument("-t","--ot",metavar="",type=str,help = "CpG_OT* file (top strand) from bismark methylation extractor")
parser.add_argument("-@","--threads",metavar="",type=int,help = "number of additional threads to use")
parser.add_argument("-p","--prefix",metavar="",type=str,help = "the prefix of the output file")
parser.add_argument("-l","--levels",metavar="",type=int,default=5,help = "number of methylation levels")
parser.add_argument("--cfTAPS", action="store_true", help="set this flag for cfTAPS sequencing data")

args = parser.parse_args()

input_path = args.input
region_path = args.region
ob_path = args.ob
ot_path = args.ot
threads = args.threads
prefix = args.prefix
n_levels = args.levels
is_cfTAPS = args.cfTAPS 

level_values = ['%g' % (i/(n_levels-1)) for i in range(n_levels)]

#check input
if 'CpG_OB' not in ob_path:
    print('error:the ob bismark file is wrong!')
    exit
if 'CpG_OT' not in ot_path:
    print('error:the ot bismark file is wrong!')
    exit

#ratio assign
def ratio_assign(meth_ratio):
	return level_values[int(meth_ratio*(n_levels-1) + 0.5)]

# Merge and sort Bismark CpG reports (OB/OT), 
call("sed '1d' -i " + ob_path,shell = True)
call("sed '1d' -i " + ot_path,shell = True)
call("cat " + ob_path + " " + ot_path + " | sort -t '\t' -k 1 | uniq > " + prefix + ".tmp.txt",shell = True)
meth_dict = {}
with open(prefix + ".tmp.txt",'r') as file:
	meth_count = 0
	total_count = 0
	prev = ""
	for line in file:
		line = line.strip().split()
		id = line[0]
		meth = line[1]
		if prev != id:
			if total_count >=3:
				ratio = meth_count / total_count
				if is_cfTAPS:
					ratio = 1 - ratio
				meth_dict[prev] = ratio_assign(ratio)
			meth_count = 0
			total_count = 0
			prev = id
		total_count += 1
		if meth == "+":
			meth_count += 1

# Check if a read belongs to a region based on overlap proportion
def reads_select(region_start,region_end,read_start,read_end,proportion):
	read_length = int(read_end) - int(read_start)
	if (read_start >= region_start) & (read_end <= region_end):
		read_inter_length = read_length
	elif (read_start > region_start) & (read_end > region_end):
		read_inter_length = region_end - read_start
	elif (region_start > read_start) & (region_end > read_end):
		read_inter_length = read_end - region_start
	elif (region_start > read_start) & (read_end > region_end):
		read_inter_length = region_end -region_start
	return read_inter_length/read_length >= proportion

#call read-level methylation
call("samtools bedcov " + region_path + " " + input_path + " | awk -v OFS='\t' '$5>0{print $1,$2,$3,$4}' > " + prefix + ".tmp.bed",shell = True)
def call_celfeer(line):
	line = line.strip().split()
	region_chr = line[0]
	region_start = int(line[1])
	region_end = int(line[2])
	region = line[3]
	bam = pysam.AlignmentFile(input_path,'rb')
	region_dict = {}
	region_dict[region] = {lv:0 for lv in level_values}
	for read in bam.fetch(contig=region_chr,start=region_start,end=region_end):
		if read.isize > 0:
			istart = read.reference_start
			iend = read.reference_start + read.isize
			select_1 =  reads_select(region_start,region_end,istart,iend,0.5)
			select_2 = read.query_name in meth_dict.keys()
			if select_1 & select_2:
				region_dict[region][meth_dict[read.query_name]] += 1
	return region_dict

#parallel computing
if __name__ == '__main__':
	with open(prefix + ".tmp.bed","r") as file:
		pool = Pool(threads)
		region_list = pool.map(call_celfeer,file)
		pool.close()
		pool.join()

#output result file
with open(prefix + ".csv",'w') as out_file:
	out_file.write("region," + ",".join(level_values) + "\n")
	for region_dict in region_list:
		region = list(region_dict.keys())[0]
		if not any(region_dict[region].values()):
			continue
		else:
			out_file.write(region + "," + ",".join(str(region_dict[region][lv]) for lv in level_values) + "\n")
out_file.close()
call("rm " + prefix + ".tmp*",shell = True)