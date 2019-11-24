from Bio.Phylo.PAML import codeml

cml = codeml.Codeml(alignment = "nonxtome_codonaligned.phy", tree = "./phylogenies/nonxtome_codon.new.treefile",
                    out_file = "results.out", working_dir = "./paml")

cml.run()
