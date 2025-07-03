import os
import sys

if len(sys.argv) < 2:
    print("Expecting input...")
    sys.exit()
elif len(sys.argv) == 2:
    input_file = sys.argv[1]
    if input_file[-4:] != ".tsv":
        print("Error: Expected a .tsv file as input.")
        sys.exit()
    output_file = input_file[:-4] + ".populated" + input_file[-4:]
elif len(sys.argv) == 3:
    input_file = sys.argv[1]
    if input_file[-4:] != ".tsv":
        print("Error: Expected a .tsv file as input.")
        sys.exit()
    output_file = sys.argv[2]
    if output_file[-4:] != ".tsv":
        print("Error: Expected a .tsv file as output.")
        sys.exit()
else:
    print("Too many inputs.")
    sys.exit()

#print(input_file)
#print(new_file_name)
with open(input_file, "r", encoding = "windows-1252") as og_file:
    with open(output_file, "w", encoding = "utf-8") as new_file:
        header = 1
        exomiser_index = 7
        for line in og_file:
            line_split = line.strip().split("\t")
            if header:
                header = 0
                for n_h in range(len(line_split)):
                    if line_split[n_h] == "Exomiser":
                        exomiser_index = n_h
                        break
                new_line = "\t".join(line_split[:exomiser_index])+"\tRANK\tID\tGENE_SYMBOL\tENTREZ_GENE_ID\tMOI\tP-VALUE\tEXOMISER_GENE_COMBINED_SCORE\tEXOMISER_GENE_PHENO_SCORE\tEXOMISER_GENE_VARIANT_SCORE\tEXOMISER_VARIANT_SCORE\tCONTRIBUTING_VARIANT\tWHITELIST_VARIANT\tFUNCTIONAL_CLASS\tHGVS\tEXOMISER_ACMG_CLASSIFICATION\tEXOMISER_ACMG_EVIDENCE\tEXOMISER_ACMG_DISEASE_ID\tEXOMISER_ACMG_DISEASE_NAME\t"+"\t".join(line_split[exomiser_index+1:])
                #print(new_line)
                new_file.write(new_line+"\n")
            else:
                exomiser_output = line_split[exomiser_index][1:-1].split("},{")
                for mode in exomiser_output:
                    new_line = "\t".join(line_split[:exomiser_index])+"\t"+"\t".join(mode.split("|"))+"\t"+"\t".join(line_split[exomiser_index+1:])
                    #print(new_line)
                    new_file.write(new_line+"\n")
