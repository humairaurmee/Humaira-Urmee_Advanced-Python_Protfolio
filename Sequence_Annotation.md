```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")

```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Computational group"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational group
    


```python
simple_seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational group', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])

```

    None. This is just an example
    


```python
simple_seq_r

```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational group', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred _quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred _quality': [40, 40, 38, 30]}
    


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816.1', description='NC_005816.1 Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.description
```




    'NC_005816.1 Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
recordGB = SeqIO.read("NC_005816.gb", "genbank")
```


```python
recordGB
```




    SeqRecord(seq=Seq(None, length=9609), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['BioProject:PRJNA224116', 'BioSample:SAMN02602970', 'Assembly:GCF_000007885.1'])




```python
recordGB.id


```




    'NC_005816.1'




```python
recordGB.name


```




    'NC_005816'




```python
recordGB.letter_annotations


```




    {}




```python
recordGB.annotations

len(recordGB.annotations)

recordGB.annotations["source"]

recordGB.dbxrefs

recordGB.features

len(recordGB.features)
```




    21




```python
len(recordGB.annotations)
```




    14




```python
recordGB.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
recordGB.dbxrefs
```




    ['BioProject:PRJNA224116',
     'BioSample:SAMN02602970',
     'Assembly:GCF_000007885.1']




```python
recordGB.features
```




    [SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(9609), strand=1), type='source', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(86), ExactPosition(1109), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(86), ExactPosition(1109), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(1105), ExactPosition(1888), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(1105), ExactPosition(1888), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(2924), ExactPosition(3119), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(2924), ExactPosition(3119), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(3784), ExactPosition(3898), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(3784), ExactPosition(3898), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(4384), ExactPosition(4780), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(4384), ExactPosition(4780), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(4814), ExactPosition(5888), strand=-1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(4814), ExactPosition(5888), strand=-1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(6115), ExactPosition(6421), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(6115), ExactPosition(6421), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(6663), ExactPosition(7602), strand=1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(6663), ExactPosition(7602), strand=1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(7788), ExactPosition(8088), strand=-1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(7788), ExactPosition(8088), strand=-1), type='CDS', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(8087), ExactPosition(8429), strand=-1), type='gene', qualifiers=...),
     SeqFeature(SimpleLocation(ExactPosition(8087), ExactPosition(8429), strand=-1), type='CDS', qualifiers=...)]




```python
len(recordGB.features)
```




    21




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9,left=8, right=9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]
    


```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.start)
```




    5




```python
int(my_location.end)
```




    9




```python
exact_location = SeqFeature.SimpleLocation(5, 9)
```


```python
exact_location.start
```




    ExactPosition(5)




```python
exact_location.end
```




    ExactPosition(9)




```python
rec1 = SeqRecord(
    Seq(
        "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
        "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
        "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
        "SSAC",
    ),
    id="gi|14150838|gb|AAK54648.1|AF376133_1",
    description="chalcone synthase [Cucumis sativus]",
)
```


```python
print(rec1)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')
    


```python
recordGB
```




    SeqRecord(seq=Seq(None, length=9609), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['BioProject:PRJNA224116', 'BioSample:SAMN02602970', 'Assembly:GCF_000007885.1'])




```python
len(recordGB)
```




    9609




```python
len(recordGB.features)
```




    21




```python
print(recordGB.features[20])
```

    type: CDS
    location: [8087:8429](-)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: inference, Value: ['COORDINATES: similar to AA sequence:RefSeq:WP_002225538.1']
        Key: locus_tag, Value: ['YP_RS22250']
        Key: note, Value: ['Derived by automated computational analysis using gene prediction method: Protein Homology.']
        Key: old_locus_tag, Value: ['pPCP10', 'YP_pPCP10']
        Key: product, Value: ['type II toxin-antitoxin system RelE/ParE family toxin']
        Key: protein_id, Value: ['WP_002225538.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MVLFSQRFDDWLNEQEDALQEKVLADLKKLQVYGPELPRPYADTVKGSRYKNMKELRVQFSGRPIRAFYAFDPIRRAIVLCAGDKSNDKRFYEKLVRIAEDEFTAHLNTLESK']
    
    


```python
print(recordGB.features[15])
```

    type: gene
    location: [6663:7602](+)
    qualifiers:
        Key: gene, Value: ['pla']
        Key: locus_tag, Value: ['YP_RS22240']
        Key: old_locus_tag, Value: ['pPCP08', 'YP_pPCP08']
    
    


```python
sub_recordGB = recordGB[4300:4800]
```


```python
len(sub_recordGB)
```




    500




```python
print(sub_recordGB.features)
```

    [SeqFeature(SimpleLocation(ExactPosition(84), ExactPosition(480), strand=1), type='gene', qualifiers=...), SeqFeature(SimpleLocation(ExactPosition(84), ExactPosition(480), strand=1), type='CDS', qualifiers=...)]
    


```python
len(sub_recordGB.features)
```




    2




```python
sub_recordGB.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(84), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_recordGB.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(84), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_recordGB.features[0])
```

    type: gene
    location: [84:480](+)
    qualifiers:
        Key: locus_tag, Value: ['YP_RS22225']
        Key: old_locus_tag, Value: ['pPCP05', 'YP_pPCP05']
    
    


```python
print(sub_recordGB.features[1])
```

    type: CDS
    location: [84:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: inference, Value: ['COORDINATES: similar to AA sequence:RefSeq:WP_001247799.1']
        Key: locus_tag, Value: ['YP_RS22225']
        Key: note, Value: ['Derived by automated computational analysis using gene prediction method: Protein Homology.']
        Key: old_locus_tag, Value: ['pPCP05', 'YP_pPCP05']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['WP_223661135.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    
    


```python
sub_recordGB.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_recordGB.dbxrefs
```




    []




```python
sub_recordGB.annotations["topology"] = "linear"
```


```python
sub_recordGB.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_recordGB.id
```




    'NC_005816.1'




```python
sub_recordGB.name
```




    'NC_005816'




```python
sub_recordGB.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_recordGB.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python

```
