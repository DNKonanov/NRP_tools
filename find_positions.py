from Bio.SeqIO import parse
from Bio.Seq import translate, reverse_complement

import warnings
warnings.filterwarnings('ignore')

def gen_ex_adomains_file(adoms_file, genome_file):
    
    
    genome = parse(genome_file, format='fasta')
    adoms = parse(adoms_file, format='fasta')

    
    adoms_ex = open('adom.ex.faa', 'w')
    for l in genome:
        genome_seq = l.seq
        break

    for domain in adoms:
        print(domain.description)
        header = '>*{}*|*'.format(domain.description)

        for i in range(3):
            loc = translate(genome_seq[i:]).find(domain.seq)
            if loc != -1:
                header += '{}:{}*+'.format(
                    loc+i, loc+i+3*len(domain.seq))

        for i in range(3):
            loc = translate(reverse_complement(genome_seq)[i:]).find(domain.seq)
            if loc != -1:
            
                header += '{}:{}*-'.format(
                        len(genome_seq) - (loc+i),
                        len(genome_seq) - (loc+i) - 3*len(domain.seq),
                    )

        adoms_ex.write(header)
        adoms_ex.write('\n')
        adoms_ex.write(str(domain.seq))
        adoms_ex.write('\n')

    adoms_ex.close()

        

    
