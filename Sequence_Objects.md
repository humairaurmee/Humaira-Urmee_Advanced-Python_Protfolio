# Unit1
## Sequence Objects
```python
from Bio.Seq import Seq
```


```python
# Create a DNA sequence
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s"% (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G
    


```python
print(len(my_seq))
```

    5
    


```python
print(my_seq[0])
```

    G
    


```python
print(my_seq[4])
```

    G
    


```python
print(my_seq[2])
```

    T
    


```python
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
len(my_seq)
```




    32




```python
my_seq.count("G")
```




    9




```python
100 * (my_seq.count("G") + my_seq.count("C"))/ len(my_seq)
```




    46.875




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
gc_fraction(my_seq)
```




    0.46875




```python
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
my_seq[2:3]
```




    Seq('T')




```python
str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    
    


```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" * 10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
protein_seq = Seq("EVRNAK")
```


```python
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table=2)
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(to_stop=True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table=2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table=2, stop_symbol="!")
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAGGTGCAATCATGTGTACTGCCTTCCTGGTTCGTGTCGCCTCATGGAAGCACAGGCCTGGAATTACTGTGCGTCATGAATACAGTAAAGGCATCGGTGATAATGCCTTATTACGTGGAGTGCGCATCGCGCGCGCATGGTGTGGAAACGATTATGAATGGCCTGCGTCTGACCATGCCTACAGTGCTGGAAGTGCAGGCACATAAGAAAGTGCCTCATCATCATCGGTGGCAGCAGCATGACACCCGCCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKRCNHVYCLPGSCRLMEAQAWNYCAS*IQ*RHR**CLITWSAHRARMVWKRL*...PA*')




```python
gene.translate(table = "Bacterial", to_stop= True)
```




    Seq('VKRCNHVYCLPGSCRLMEAQAWNYCAS')




```python
gene.translate(table = "Bacterial", cds= False)
```




    Seq('VKRCNHVYCLPGSCRLMEAQAWNYCAS*IQ*RHR**CLITWSAHRARMVWKRL*...PA*')




```python
from Bio.Data import CodonTable
```


```python
#print(standard_table)
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAACCTGAATGTGAGGTGAGGTCAAGGATAGT"}, length=159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CTGAATGTGA')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAACCTGAATGTGAGGTGAGGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length=10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTATAGCGGCCGTAAAGGTGCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTATAGCGGCCGTAAAGGTGCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTATAGCGGCCGTAAAGGTGCCGA')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCGTGGAAATGCCGGCGATATGCTACCG')




```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCGTGGAAATGCCGGCGATATGCTACCG')




```python

```
