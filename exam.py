import subprocess
import os
import sys
from Bio import SeqIO

def is_fasta(filename):
    with open(filename,"r") as handle:
        fasta = SeqIO.parse(handle,"fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def run_cmd(cmd,show=True):
	#use to run command without stdout
	#inout command line

	#if need to print command to screen , show=True
	if show:
		print(cmd)
	process=subprocess.call(cmd, shell=True)

def build_db(input,type,output):
    #build database
    mkdb_cmd="makeblastdb -in {} -dbtype {}  -parse_seqids -out {} -logfile makeblastdb.log".format(input,type,output)
    run_cmd(mkdb_cmd,show_command)

def run_software(software,db,query_seq):
    #select software and run it
    if software=="blastp":
        blastp_cmd="blastp -db {} -query {} -outfmt 7 > blastp.output".format(db,query_seq)
        run_cmd(blastp_cmd,True)
        print("Run blastp finished , output file is blastp.output")
    if software=="blastn":
        blastn_cmd="blastn -db {} -query {} -outfmt 7 > blastn.output".format(db,query_seq)
        run_cmd(blastn_cmd,True)
        print("Run blastn finished , output file is blastn.output")
    if software=="blastx":
        blastx_cmd="blastx -db {} -query {} -outfmt 7 > blastx.output".format(db,query_seq)
        run_cmd(blastx_cmd,True)
        print("Run blastx finished , output file is blastx.output")
    if software=="tblastn":
        tblastn_cmd="tblastn -db {} -query {} -outfmt 7 > tblastn.output".format(db,query_seq)
        run_cmd(tblastn_cmd,True)
        print("Run tblastn finished , output file is tblastn.output")
    if software=="megablast":
        blastn_cmd="blastn -task megablast -db {} -query {} -outfmt 7 > megablast.output".format(db,query_seq)
        run_cmd(blastn_cmd,True)
        print("Run megablast finished , output file is megablast.output")
    if software=="tblastx":
        tblastx_cmd="tblastx -db {} -query {} -outfmt 7 > tblastx.output".format(db,query_seq)
        run_cmd(tblastx_cmd,True)
        print("Run tblastx finished , output file is tblastx.output")
    if software=="psi-blast":
        num_iterations=input("num_iterations[integer]: \n")
        max_target_seqs=input("max_target_seqs[integer]: \n")
        if num_iterations.isdigit() and max_target_seqs.isdigit():
            psiblast_cmd="psiblast -db {} -query {}  -num_iterations {} -seg yes -outfmt 7 -max_target_seqs {} > psiblast.output".format(db,query_seq,num_iterations,max_target_seqs)
        else:
            print("Can not recognize the parameter, using num_iterations=3 and max_target_seqs=100")
            psiblast_cmd="psiblast -db {} -query {}  -num_iterations 3 -seg yes -outfmt 7 -max_target_seqs 100 > psiblast.output".format(db,query_seq)
        run_cmd(psiblast_cmd,True)
        print("Run psi-blast finished , output file is psiblast.output.output")

#input of query type and path
while True:
    query_type=input("the type of query [nucl|prot]:\n")
    if query_type in ["nucl","prot"]:
        break
    else:
        print("Please input \"nucl\" or \"prot\"")
while True:
    query_seq=input("the path of your query file [fasta formate]:\n")
    if os.path.exists(query_seq):
        if is_fasta(query_seq):
            break
        else:
            print("the type of your sequence is not fasta formate")
            sys.exit()
    else:
        print("There is NOT query file in path:{} \n please input correct path\n if you want to quit using Ctrl+C ".format(query_seq))

#input of database type and path
while True:
    database_type=input("the type of query [nucl|prot]:\n")
    if database_type in ["nucl","prot"]:
        break
    else:
        print("Please input \"nucl\" or \"prot\"")

#choose database or create new database
database_exist=input("do you want to creat a new blast data base? [Y|N]\n")
if database_exist=="N":
        while True:
            #differnet database have differnet suffix prot database is .pdb nucl is .ndb
            database_db=input("the path of your database:\n")
            if database_type=="nucl":
                test=database_db+".ndb"
            else:
                test=database_db+".pdb"
            if os.path.exists(test):
                break
            else:
                print("There is NOT Corresponding database in path:{}\n please input correct path\n if you want to quit using Ctrl+C ".format(database_db))
elif database_exist=="Y":
    while True:
        database_seq=input("the path of your database sequence[fasta formate]:\n")
        if os.path.exists(database_seq):
            if is_fasta(database_seq):
                try:
                    os.makedirs("database")
                except:
                    pass
                build_db(database_seq,database_type,"./database/{}_database".format(database_type))
                database_db="./database/{}_database".format(database_type)
                break
            else:
                print("the type of your sequence is not fasta\n please check formate \n if you want to quit using Ctrl+C")
        else:
            print("There is NOT database sequence in path:{} \n please check formate \n if you want to quit using Ctrl+C".format(database_seq))
else:
    print("Please input \"Y\" or \"N\"")
    sys.exit()

#show the corresponding software for query type and database type
query_database=(query_type+"_"+database_type)
software_dict={"nucl_prot":["blastx"],"nucl_nucl":["blastn","megablast","tblastx"],"prot_prot":["blastp","psi-blast"],"prot_nucl":["tblastn"]}

print("Your query type is {}， and the database type is {} .\n".format(query_type,database_type))
software_can_use=software_dict.get(query_database)
print("these software can use:{}".format(software_can_use))

while True:
    which_command=input("select the software you want to use:\n")
    if which_command in software_can_use:
        break
    else:
        print("Your input not in list, please input again")

run_software(which_command,database_db,query_seq)

