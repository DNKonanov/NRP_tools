// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            National Center for Biotechnology Information (NCBI)
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act. It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted. This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government do not place any restriction on its use or reproduction.
//  We would, however, appreciate having the NCBI and the author cited in
//  any work or product based on this material.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
// ===========================================================================
//
// File Name:  xtract.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"bufio"
	"bytes"
	"container/heap"
	"fmt"
	"html"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
	"unicode"
)

// VERSION AND HELP MESSAGE TEXT

const xtract_version = "4.20"

const xtract_help = `
Exploration Argument Hierarchy

  -pattern         Name of record within set
  -group             Use of different argument
  -block               names allows command-line
  -subset                control of nested looping

Conditional Execution

  -position        Must be at given location in list
  -match           Element [@attribute] [:value] required
  -avoid           Skip if element matches
  -and             All tests must pass
  -or              Any passing test suffices

Format Customization

  -ret             Override line break between patterns
  -tab             Replace tab character between fields
  -sep             Separator between group members
  -pfx             Prefix to print before group
  -sfx             Suffix to print after group
  -lbl             Insert arbitrary text

Item Selection

  -element         Print all items that match tag name
  -first           Only print value of first item
  -last            Only print value of last item
  -NAME            Record value in named variable

-element Constructs

  Tag              Caption
  Group            Initials,LastName
  Parent/Child     MedlineCitation/PMID
  Attribute        DescriptorName@MajorTopicYN
  Object Count     "#Author"
  Item Length      "%Title"
  Variable         "&NAME"

Exploration Constructs

  Object           DateCreated
  Parent/Child     Book/AuthorList
  Heterogeneous    PubmedArticleSet/*
  Recursive        */Taxon

Command Generator

  -insd            Generate INSDSeq extraction commands

-insd Argument Order

  Descriptors      INSDSeq_sequence INSDSeq_definition INSDSeq_division
  Flags            complete or partial [optional]
  Feature(s)       CDS,mRNA
  Qualifiers       INSDFeature_key "#INSDInterval" gene product

XML Formatting

  -format          Repair XML format and indentation
  -outline         Display outline of XML structure
  -synopsis        Display count of unique XML paths

Documentation

  -help            Print this document
  -extras          Additional constraint and selection commands
  -examples        Examples of EDirect and xtract usage
  -scripts         Advanced EDirect and xtract automation examples
  -version         Print version number

Examples

  -pattern DocumentSummary -element Id -first Name Title

  -pattern "PubmedArticleSet/*" -block Author -sep " " -element Initials,LastName

  -pattern PubmedArticle -block MeshHeading -match "@MajorTopicYN:Y" -sep " / " -element DescriptorName,QualifierName

  -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop

  -pattern Taxon -block "*/Taxon" -avoid "Rank:no rank" -tab "\n" -element Rank,ScientificName

  -insd CDS gene product protein_id translation

  -insd complete mat_peptide "%peptide" product peptide
`

const xtract_extras = `
XML Processing

  -cleanup        Fix non-ASCII spaces
  -compress       Compress runs of spaces

String -match Constraints (Case-Insensitive)

  -equals         String must match exactly
  -contains       Substring must be present
  -starts-with    Substring must be at beginning
  -ends-with      Substring must be at end

Numeric -match Constraints (Integer Values)

  -gt             Greater than
  -ge             Greater than or equal to
  -lt             Less than
  -le             Less than or equal to
  -eq             Equal to
  -ne             Not equal to

Format Customization

  -clr            Clear queued tab separator
  -pfc            Preface combines -clr and -pfx
  -rst            Reset -sep, -pfx, and -sfx
  -head           Print text before all results
  -tail           Print text after all results

Item Selection

  -encode         URL-encode <, >, &, ", and ' characters

Numeric Selection

  -sum            Sum
  -min            Minimum
  -max            Maximum

-element Constructs

  Item Depth      "^PMID"
  Parent Index    "+"
  XML Subtree     "*"

XML Validation

  -verify         Report XML data integrity problems

Examples

  -match "%Title" -le 70 -and "#Author" -lt 6

  -match "^PMID" -eq 3

  -match "DateCreated/Year" -ne "DateRevised/Year"

  -match CommonName -contains mouse

  -match "&ABST" -starts-with "Transposable elements"

  -block INSDReference -position 2

  -min ChrStart,ChrStop

  -max ExonCount
`

const xtract_internal = `
Legacy Parsing Method

  -legacy    Parser reads input token stream

Performance Tuning

  -proc      Maximum number of CPU processors
  -serv      Number of concurrent parser instances
  -chan      Depth of communication channels
  -heap      Size of unshuffling heap

Debugging

  -timer     Report processing duration and rate
  -debug     Display run-time parameters
  -empty     Flag records with no output
  -index     Print record index numbers

Internal Performance Tests

  -chunk     ReadBlocks
  -split     ReadBlocks -> SplitPattern
  -token     ReadBlocks -> SplitPattern -> StreamTokens
  -parse     ReadBlocks -> SplitPattern -> StreamTokens -> ParseXML
  -trial     ReadBlocks -> SplitPattern -> StreamTokens -> ParseXML -> ProcessQuery

Documentation

  -keys      Keyboard navigation shortcuts
  -unix      Common Unix commands
`

const xtract_examples = `
Publications

  efetch -db pubmed -id 6271474,5685784,4882854,6243420 -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID "#Author" \
    -block Author -position first -sep " " -element Initials,LastName \
    -block Article -element ArticleTitle

  6271474    5    MJ Casadaban     Tn3: transposition and control.
  5685784    2    RK Mortimer      Suppressors and suppressible mutations in yeast.
  4882854    2    ED Garber        Proteins and enzymes as taxonomic tools.
  6243420    1    NR Cozzarelli    DNA gyrase and the supercoiling of DNA.

Formatted Authors

  efetch -db pubmed -id 1413997,6301692,781293 -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID \
    -block DateCreated -sep "-" -element Year,Month,Day \
    -block Author -sep " " -tab "" \
      -element "&COM" Initials,LastName -COM "(, )"

  1413997    1992-11-25    RK Mortimer, CR Contopoulou, JS King
  6301692    1983-06-17    MA Krasnow, NR Cozzarelli
  781293     1976-10-02    MJ Casadaban

Medical Subject Headings

  efetch -db pubmed -id 6092233,2539356,1937004 -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID \
    -block MeshHeading \
      -subset DescriptorName -pfc "\n" -sep "|" -element @MajorTopicYN,DescriptorName \
      -subset QualifierName -pfc " / " -sep "|" -element @MajorTopicYN,QualifierName |
  sed -e 's/N|//g' -e 's/Y|/*/g'

  6092233
  Base Sequence
  DNA Restriction Enzymes
  DNA, Fungal / genetics / *isolation & purification
  *Genes, Fungal
  ...

Peptide Sequences

  esearch -db protein -query "conotoxin AND mat_peptide [FKEY]" |
  efetch -format gpc |
  xtract -insd complete mat_peptide "%peptide" product peptide |
  grep -i conotoxin | sort -t $'\t' -u -k 2,2n | head -n 8

  ADB43131.1    15    conotoxin Cal 1b     LCCKRHHGCHPCGRT
  AIC77099.1    16    conotoxin Im1.2      GCCSHPACNVNNPHIC
  AIC77105.1    17    conotoxin Lt1.4      GCCSHPACDVNNPDICG
  AIC77103.1    18    conotoxin Lt1.2      PRCCSNPACNANHAEICG
  AIC77083.1    20    conotoxin Bt14.6     KDCTYCMHSSCSMMYEKCRP
  AIC77085.1    21    conotoxin Bt14.8     NECDNCMRSFCSMIYEKCRLK
  AIC77093.1    22    conotoxin Bt14.16    GDCKPCMHPDCRFNPGRCRPRE
  AIC77154.1    23    conotoxin Bt14.19    VREKDCPPHPVPGMHKCVCLKTC

Chromosome Locations

  esearch -db gene -query "calmodulin [PFN] AND mammalia [ORGN]" |
  efetch -format docsum |
  xtract -pattern DocumentSummary -MAP "(-)" -MAP MapLocation \
    -element Id Name "&MAP" ScientificName

  801       CALM1    14q32.11         Homo sapiens
  808       CALM3    19q13.2-q13.3    Homo sapiens
  805       CALM2    2p21             Homo sapiens
  24242     Calm1    6q31-q32         Rattus norvegicus
  12313     Calm1    12 E             Mus musculus
  326597    CALM     -                Bos taurus
  50663     Calm2    6q11-q12         Rattus norvegicus
  24244     Calm3    1q22             Rattus norvegicus
  12315     Calm3    7 9.15 cM        Mus musculus
  12314     Calm2    17 E4            Mus musculus
  617095    CALM1    -                Bos taurus
  396838    CALM3    6                Sus scrofa
  ...

Gene Regions

  esearch -db gene -query "DDT [GENE] AND mouse [ORGN]" |
  efetch -format docsum |
  xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop |
  xargs -n 3 sh -c 'efetch -db nuccore -format gb \
    -id "$0" -chr_start "$1" -chr_stop "$2"'

  LOCUS       NC_000076               2142 bp    DNA     linear   CON 09-FEB-2015
  DEFINITION  Mus musculus strain C57BL/6J chromosome 10, GRCm38.p3 C57BL/6J.
  ACCESSION   NC_000076 REGION: complement(75771233..75773374) GPC_000000783
  VERSION     NC_000076.6  GI:372099100
  ...
  FEATURES             Location/Qualifiers
       source          1..2142
                       /organism="Mus musculus"
                       /mol_type="genomic DNA"
                       /strain="C57BL/6J"
                       /db_xref="taxon:10090"
                       /chromosome="10"
       gene            1..2142
                       /gene="Ddt"
       mRNA            join(1..159,462..637,1869..2142)
                       /gene="Ddt"
                       /product="D-dopachrome tautomerase"
                       /transcript_id="NM_010027.1"
       CDS             join(52..159,462..637,1869..1941)
                       /gene="Ddt"
                       /codon_start=1
                       /product="D-dopachrome decarboxylase"
                       /protein_id="NP_034157.1"
                       /translation="MPFVELETNLPASRIPAGLENRLCAATATILDKPEDRVSVTIRP
                       GMTLLMNKSTEPCAHLLVSSIGVVGTAEQNRTHSASFFKFLTEELSLDQDRIVIRFFP
                       ...

Taxonomic Names

  esearch -db taxonomy -query "txid10090 [SBTR] OR camel [COMN]" |
  efetch -format docsum |
  xtract -pattern DocumentSummary -match CommonName \
    -element Id ScientificName CommonName

  57486    Mus musculus molossinus    Japanese wild mouse
  39442    Mus musculus musculus      eastern European house mouse
  35531    Mus musculus bactrianus    southwestern Asian house mouse
  10092    Mus musculus domesticus    western European house mouse
  10091    Mus musculus castaneus     southeastern Asian house mouse
  10090    Mus musculus               house mouse
  9838     Camelus dromedarius        Arabian camel
  9837     Camelus bactrianus         Bactrian camel

Structural Similarity

  esearch -db structure -query "crotalus [ORGN] AND phospholipase A2" |
  elink -related | efilter -query "archaea [ORGN]" | efetch -format docsum |
  xtract -pattern DocumentSummary -match "PdbClass:Hydrolase" \
    -element PdbAcc PdbDescr

  3VV2    Crystal Structure Of Complex Form Between S324a-subtilisin And Mutant Tkpro
  3VHQ    Crystal Structure Of The Ca6 Site Mutant Of Pro-Sa-Subtilisin
  2ZWP    Crystal Structure Of Ca3 Site Mutant Of Pro-S324a
  2ZWO    Crystal Structure Of Ca2 Site Mutant Of Pro-S324a
  ...

Indexed Fields

  einfo -db pubmed |
  xtract -pattern Field -match "IsDate:Y" -and "IsHidden:N" \
    -pfx "[" -sep "]\t" -element Name,FullName |
  sort -t $'\t' -k 2f

  [CDAT]    Date - Completion
  [CRDT]    Date - Create
  [EDAT]    Date - Entrez
  [MHDA]    Date - MeSH
  [MDAT]    Date - Modification
  [PDAT]    Date - Publication

Record Counts

  echo "diphtheria measles pertussis polio tuberculosis" |
  xargs -n 1 sh -c 'esearch -db pubmed -query "$0 [MESH]" |
  efilter -days 365 -datetype PDAT |
  xtract -pattern ENTREZ_DIRECT -lbl "$0" -element Count'

  diphtheria      18
  measles         166
  pertussis       98
  polio           75
  tuberculosis    1386
`

const xtract_scripts = `
Gene Products

  for sym in HBB DMD TTN ATP7B HFE BRCA2 CFTR PAH PRNP RAG1
  do
    esearch -db gene -query "$sym [GENE] AND human [ORGN]" |
    efilter -query "alive [PROP]" | efetch -format docsum |
    xtract -pattern GenomicInfoType \
      -element ChrAccVer ChrStart ChrStop |
    while read acc str stp
    do
      efetch -db nuccore -format gbc \
        -id "$acc" -chr_start "$str" -chr_stop "$stp" |
      xtract -insd CDS,mRNA INSDFeature_key "#INSDInterval" \
        gene "%transcription" "%translation" \
        product transcription translation |
      grep -i $'\t'"$sym"$'\t'
    done
  done

  NC_000011.10    mRNA    3     HBB    626      hemoglobin, beta                     ACATTTGCTT...
  NC_000011.10    CDS     3     HBB    147      hemoglobin subunit beta              MVHLTPEEKS...
  NC_000023.11    mRNA    78    DMD    13805    dystrophin, transcript variant X2    AGGAAGATGA...
  NC_000023.11    mRNA    77    DMD    13794    dystrophin, transcript variant X6    ACTTTCCCCC...
  NC_000023.11    mRNA    77    DMD    13800    dystrophin, transcript variant X5    ACTTTCCCCC...
  NC_000023.11    mRNA    77    DMD    13785    dystrophin, transcript variant X7    ACTTTCCCCC...
  NC_000023.11    mRNA    74    DMD    13593    dystrophin, transcript variant X8    ACTTTCCCCC...
  NC_000023.11    mRNA    75    DMD    13625    dystrophin, transcript variant X9    ACTTTCCCCC...
  ...

Genome Range

  BetweenTwoGenes() {
    awk -F '\t' -v 'OFS=\t' \
      "/^$1\t/{a++}/^$2\t/{a++}a>0{print}a>1{exit}"
  }

  esearch -db gene -query "Homo sapiens [ORGN] AND 13 [CHR]" |
  efilter -query "alive [PROP]" | efetch -format docsum |
  xtract -pattern DocumentSummary -NAME Name -DESC Description \
    -block GenomicInfoType -match "ChrLoc:13" \
      -element "&NAME" "&DESC" -tab "\n" -min ChrStart,ChrStop |
  sort -t $'\t' -k 3,3n | cut -f 1-2 | BetweenTwoGenes NEK5 DHRS12

  DHRS12       dehydrogenase/reductase (SDR family) member 12
  LINC00282    long intergenic non-protein coding RNA 282
  CCDC70       coiled-coil domain containing 70
  ATP7B        ATPase, Cu++ transporting, beta polypeptide
  CTAGE3P      CTAGE family member 3, pseudogene
  FABP5P2      fatty acid binding protein 5 pseudogene 2
  ALG11        ALG11, alpha-1,2-mannosyltransferase
  UTP14C       UTP14, U3 small nucleolar ribonucleoprotein, homolog C (yeast)
  NEK5         NIMA-related kinase 5

Amino Acid Substitutions

  ApplySNPs() {
    seq=""
    last=""

    while read rsid accn pos res
    do
      if [ "$accn" != "$last" ]
      then
        insd=$(efetch -db protein -id "$accn" -format gbc < /dev/null)
        seq=$(echo $insd | xtract -pattern INSDSeq -element INSDSeq_sequence)
        last=$accn
      fi

      pos=$((pos+1))
      pfx=""
      sfx=""
      echo ">rs$rsid [$accn $res@$pos]"
      if [ $pos -gt 1 ]
      then
        pfx=$(echo ${seq:0:$pos-1})
      fi
      if [ $pos -lt ${#seq} ]
      then
        sfx=$(echo ${seq:$pos})
      fi
      echo "$pfx$res$sfx" | fold -w 50
    done
  }

  esearch -db gene -query "OPN1MW [GENE] AND human [ORGN]" |
  elink -target snp | efetch -format xml |
  xtract -pattern Rs -RSID Rs@rsId \
    -block FxnSet -match @fxnClass:missense \
      -sep "." -element "&RSID" @protAcc,@protVer @aaPosition \
      -tab "\n" -element @residue |
  sort -t $'\t' -k 2,2 -k 3,3n -k 4,4 | uniq | ApplySNPs

  >rs104894915 [NP_000504.1 K@94]
  maqqwslqrlagrhpqdsyedstqssiftytnsnstrgpfegpnyhiapr
  wvyhltsvwmifvviasvftnglvlaatmkfkklrhplnwilvKlavadl
  aetviastisvvnqvygyfvlghpmcvlegytvslcgitglwslaiiswe
  ...

Amino Acid Composition

  #!/usr/bin/env bash

  abbrev=( Ala Asx Cys Asp Glu Phe Gly His Ile \
           Xle Lys Leu Met Asn Pyl Pro Gln Arg \
           Ser Thr Sec Val Trp Xxx Tyr Glx )

  AminoAcidComp() {
    local count
    while read num lttr
    do
      idx=$(printf %i "'$lttr'")
      ofs=$((idx-97))
      count[$ofs]="$num"
    done <<< "$1"
    for i in {0..25}
    do
      echo -e "${abbrev[$i]}\t${count[$i]-0}"
    done |
    sort
  }

  AminoAcidJoin() {
    result=""
    while read acc seq gene
    do
      comp="$(echo "$seq" | fold -w 1 | sort-uniq-count)"
      current=$(AminoAcidComp "$comp")
      current=$(echo -e "GENE\t$gene\n$current")
      if [ -n "$result" ]
      then
        result=$(join -t $'\t' <(echo "$result") <(echo "$current"))
      else
        result=$current
      fi
    done
    echo "$result" |
    grep -e "GENE" -e "[1-9]"
  }

  ids="NP_001172026,NP_000509,NP_004001,NP_001243779"
  efetch -db protein -id "$ids" -format gpc |
  xtract -insd INSDSeq_sequence CDS gene | grep '.' | AminoAcidJoin

  GENE    INS    HBB    DMD    TTN
  Ala     10     15     210    2084
  Arg     5      3      193    1640
  Asn     3      6      153    1111
  Asp     2      7      185    1720
  Cys     6      2      35     513
  Gln     7      3      301    942
  Glu     8      8      379    3193
  Gly     12     13     104    2066
  His     2      9      84     478
  Ile     2      0      165    2062
  Leu     20     18     438    2117
  Lys     2      11     282    2943
  Met     2      2      79     398
  Phe     3      8      77     908
  Pro     6      7      130    2517
  Ser     5      5      239    2463
  Thr     3      7      194    2546
  Trp     2      2      67     466
  Tyr     4      3      61     999
  Val     6      18     186    3184

Phrase Searching

  entrez-phrase-search -db pubmed -field WORD \
    selective serotonin reuptake inhibitor + monoamine oxidase inhibitor |
  efetch -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID \
    -block Keyword -pfc "\n  " -element Keyword

  24657329
    Antidepressant
    Organic cation transporter 2
    Piperine
    Uptake 2
  24280122
    5-HIAA
    5-HT
    5-HTP
    5-hydroxyindoleacetic acid
    5-hydroxytryptophan
    ...
`

const keyboard_shortcuts = `
Command History

  Ctrl-n     Next command
  Ctrl-p     Previous command

Move Cursor Forward

  Ctrl-e     To end of line
  Ctrl-f     By one character
  Esc-f      By one argument

Move Cursor Backward

  Ctrl-a     To beginning of line
  Ctrl-b     By one character
  Esc-b      By one argument

Delete

  Del        Previous character
  Ctrl-d     Next character
  Ctrl-k     To end of line
  Ctrl-u     Entire line
  Ctrl-w     Previous word
  Esc-Del    Previous argument
  Esc-d      Next argument

Autocomplete

  Tab        Completes directory or file names

Program Control

  Ctrl-c     Quit running program
  ^x^y       Run last command replacing x with y
`

const unix_commands = `
Process by Contents

 sort      Sorts lines of text

  -f       Ignore case
  -n       Numeric comparison
  -r       Reverse result order

  -k       Field key (start,stop or first)
  -u       Unique lines with identical keys

  -b       Ignore leading blanks
  -s       Stable sort
  -t       Specify field separator

 uniq      Removes repeated lines

  -c       Count occurrences
  -i       Ignore case

  -f       Ignore first n fields
  -s       Ignore first n characters

  -d       Only output repeated lines
  -u       Only output non-repeated lines

 grep      Matches patterns using regular expressions

  -i       Ignore case
  -v       Invert search
  -w       Search expression as a word
  -x       Search expression as whole line

  -e       Specify individual pattern

  -c       Only count number of matches
  -n       Print line numbers

Regular Expressions

 Characters

  .        Any single character (except newline)
  \w       Alphabetic [A-Za-z], numeric [0-9], or underscore (_)
  \s       Whitespace (space or tab)
  \        Escapes special characters
  []       Matches any enclosed characters

 Positions

  ^        Beginning of line
  $        End of line
  \b       Word boundary

 Repeat Matches

  ?        0 or 1
  *        0 or more
  +        1 or more
  {n}      Exactly n

 Escape Sequences

  \n       Line break
  \t       Tab character

Modify Contents

 sed       Replaces text strings

  -e       Specify individual expression

 tr        Translates characters

  -d       Delete character

 rev       Reverses characters on line

Format Contents

 column    Aligns columns by content width

  -s       Specify field separator
  -t       Create table

 expand    Aligns columns to specified positions

  -t       Tab positions

 fold      Wraps lines at a specific width

  -w       Line width

Filter by Position

 cut       Removes parts of lines

  -c       Characters to keep
  -f       Fields to keep
  -d       Specify field separator
  -s       Suppress lines with no delimiters

 head      Prints first lines

  -n       Number of lines

 tail      Prints last lines

  -n       Number of lines

Miscellaneous

 wc        Counts words, lines, or characters

  -c       Characters
  -l       Lines
  -w       Words

 xargs     Constructs arguments

  -n       Number of words per batch

Directory and File Navigation

 cd        Changes directory

  /        Root
  ~        Home
  .        Current
  ..       Parent
  -        Previous

 ls        Lists file names

  -1       One entry per line
  -a       Show files beginning with dot (.)
  -l       List in long format
  -R       Recursively explore subdirectories
  -S       Sort files by size

 pwd       Prints working directory path
`

// DATA OBJECTS

type Node struct {
	Name       string
	Parent     string
	Contents   string
	Attributes string
	Attribs    []string
	Children   *Node
	Next       *Node
}

type Operation struct {
	Type  int
	Value string
}

type Block struct {
	Visit      string
	Position   string
	Tests      []string
	Actions    []string
	Conditions []*Operation
	Commands   []*Operation
	Subtasks   []*Block
}

// UTILITIES

func IsNotJustWhitespace(str string) bool {

	for _, ch := range str {
		if ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n' && ch != '\f' {
			return true
		}
	}

	return false
}

func HasAmpOrNotASCII(str string) bool {

	for _, ch := range str {
		if ch == '&' || ch > 127 {
			return true
		}
	}

	return false
}

func IsAllCapsOrDigits(str string) bool {

	for _, rune := range str {
		if !unicode.IsUpper(rune) && !unicode.IsDigit(rune) {
			return false
		}
	}

	return true
}

func CompressRunsOfSpaces(str string) string {

	whiteSpace := false
	var buffer bytes.Buffer

	for _, rune := range str {
		if unicode.IsSpace(rune) {
			if !whiteSpace {
				buffer.WriteRune(' ')
			}
			whiteSpace = true
		} else {
			buffer.WriteRune(rune)
			whiteSpace = false
		}
	}

	return buffer.String()
}

func HasFlankingSpace(str string) bool {

	if str == "" {
		return false
	}

	if str[0] == ' ' {
		return true
	}

	strlen := len(str)
	if str[strlen-1] == ' ' {
		return true
	}

	return false
}

func HasBadSpace(str string) bool {

	for _, rune := range str {
		if unicode.IsSpace(rune) && rune != ' ' {
			return true
		}
	}

	return false
}

func CleanupBadSpaces(str string) string {

	var buffer bytes.Buffer

	for _, rune := range str {
		if unicode.IsSpace(rune) {
			buffer.WriteRune(' ')
		} else {
			buffer.WriteRune(rune)
		}
	}

	return buffer.String()
}

const (
	_ = iota
	LEFT
	RIGHT
)

func SplitInTwoAt(str, chr string, side int) (string, string) {

	slash := strings.SplitN(str, chr, 2)
	if len(slash) > 1 {
		return slash[0], slash[1]
	}

	if side == LEFT {
		return str, ""
	}

	return "", str
}

const (
	UNSET = iota
	ELEMENT
	FIRST
	LAST
	ENCODE
	PFX
	SFX
	SEP
	TAB
	RET
	LBL
	CLR
	PFC
	RST
	MATCH
	AVOID
	AND
	OR
	EQUALS
	CONTAINS
	STARTSWITH
	ENDSWITH
	GT
	GE
	LT
	LE
	EQ
	NE
	SUM
	MIN
	MAX
	VARIABLE
	VALUE
	STAR
	COUNT
	LENGTH
	DEPTH
	UNRECOGNIZED
)

func ParseFlag(str string) int {

	switch str {
	case "-element":
		return ELEMENT
	case "-first":
		return FIRST
	case "-last":
		return LAST
	case "-encode":
		return ENCODE
	case "-pfx":
		return PFX
	case "-sfx":
		return SFX
	case "-sep":
		return SEP
	case "-tab":
		return TAB
	case "-ret":
		return RET
	case "-lbl":
		return LBL
	case "-clr":
		return CLR
	case "-pfc":
		return PFC
	case "-rst":
		return RST
	case "-match":
		return MATCH
	case "-avoid":
		return AVOID
	case "-and":
		return AND
	case "-or":
		return OR
	case "-equals":
		return EQUALS
	case "-contains":
		return CONTAINS
	case "-starts-with":
		return STARTSWITH
	case "-ends-with":
		return ENDSWITH
	case "-gt":
		return GT
	case "-ge":
		return GE
	case "-lt":
		return LT
	case "-le":
		return LE
	case "-eq":
		return EQ
	case "-ne":
		return NE
	case "-sum":
		return SUM
	case "-min":
		return MIN
	case "-max":
		return MAX
	default:
		if len(str) > 1 && str[0] == '-' && IsAllCapsOrDigits(str[1:]) {
			return VARIABLE
		}
	}

	if len(str) > 0 && str[0] == '-' {
		return UNRECOGNIZED
	}

	return UNSET
}

func ConvertSlash(str string) string {

	if str == "" {
		return str
	}

	len := len(str)
	res := make([]byte, len+1, len+1)

	isSlash := false
	idx := 0
	for _, rune := range str {
		if isSlash {
			switch rune {
			case 'n':
				// line feed
				res[idx] = '\n'
			case 'r':
				// carriage return
				res[idx] = '\r'
			case 't':
				// horizontal tab
				res[idx] = '\t'
			case 'f':
				// form feed
				res[idx] = '\f'
			case 'a':
				// audible bell from terminal (undocumented)
				res[idx] = '\x07'
			default:
				res[idx] = byte(rune)
			}
			idx++
			isSlash = false
		} else if rune == '\\' {
			isSlash = true
		} else {
			res[idx] = byte(rune)
			idx++
		}
	}

	res = res[0:idx]

	return string(res)
}

// CONVERT XML INPUT DATA TO STREAM OF BLOCKS

type XmlReader struct {
	Reader     *bufio.Reader
	Buffer     []byte
	Remainder  string
	Closed     bool
	Docompress bool
	Docleanup  bool
}

func NewXmlReader(in *bufio.Reader, doCompress, doCleanup bool) *XmlReader {

	if in == nil {
		return nil
	}

	rdr := &XmlReader{Reader: in, Docompress: doCompress, Docleanup: doCleanup}

	// 65536 appears to be the maximum number of characters read by bufio.Reader
	const XMLBUFSIZE = 65536

	rdr.Buffer = make([]byte, XMLBUFSIZE)

	return rdr
}

// read buffer, concatenate if necessary to place long element content into a single string
func NextBlock(rdr *XmlReader) string {

	if rdr == nil {
		return ""
	}

	// read one buffer, trim at last > and retain remainder for next call, signal if no > character
	nextBuffer := func() (string, bool, bool) {

		if rdr.Closed {
			return "", false, true
		}

		n, err := rdr.Reader.Read(rdr.Buffer)
		if err != nil {
			// end of file
			rdr.Closed = true
			// send final remainder
			return rdr.Remainder, false, false
		}

		bufr := string(rdr.Buffer)
		// prepend previous remainder
		line := rdr.Remainder + bufr[0:n]
		rdr.Remainder = ""

		// look for last > character
		pos := -1
		for pos = len(line) - 1; pos >= 0; pos-- {
			if line[pos] == '>' {
				break
			}
		}
		if pos > -1 {
			// trim back to last > character, save remainder for next buffer
			pos++
			rdr.Remainder = line[pos:]
			line = line[0:pos]
			return line, false, false
		}

		// no > found, signal need to continue reading long content
		return line, true, false
	}

	// read next buffer
	line, cont, closed := nextBuffer()

	if closed {
		return ""
	}

	// if buffer does not end with > character
	if cont {
		var buff bytes.Buffer

		// keep reading long content blocks
		for {
			if line != "" {
				buff.WriteString(line)
			}
			if !cont {
				break
			}
			line, cont, closed = nextBuffer()
			if closed {
				break
			}
		}

		// concatenate blocks
		line = buff.String()
	}

	// trim spaces around final block
	if HasFlankingSpace(line) {
		line = strings.TrimSpace(line)
	}
	// optionally compress/cleanup tags/attributes and contents
	if rdr.Docompress {
		line = CompressRunsOfSpaces(line)
	}
	if rdr.Docleanup {
		if HasBadSpace(line) {
			line = CleanupBadSpaces(line)
		}
	}

	return line
}

// CONVERT XML BLOCKS TO STREAM OF TOKENS

const (
	_ = iota
	IS_START_TAG
	IS_CONTENT
	IS_END_TAG
	IS_SELF_TAG
	IS_NOT_A_TAG
	IS_CLOSED
)

type Tokenizer struct {
	Reader   *XmlReader
	Text     string
	TxtLen   int
	Idx      int
	Start    int
	Line     int
	EndToken string
}

func NewTokenizer(in *XmlReader, str string) *Tokenizer {

	if in == nil && str == "" {
		return nil
	}

	rdr := &Tokenizer{Reader: in, Text: str, TxtLen: len(str), Line: 1}

	return rdr
}

// return class, element name, attributes, content text
func NextToken(rdr *Tokenizer) (int, string, string, string) {

	if rdr == nil {
		return 0, "", "", ""
	}

	// function returns true for contents, false for element
	nextUnit := func() (string, bool) {

		for {
			if rdr.Text == "" && rdr.Reader != nil {
				// get next block
				rdr.Text = NextBlock(rdr.Reader)
				rdr.TxtLen = len(rdr.Text)
				rdr.Idx = 0
				rdr.Start = 0
			}

			if rdr.Text == "" {
				// if closed, break out of loop
				break
			}

			// copy current tokenizer values into local variables for speed
			text := rdr.Text[:]
			txtlen := rdr.TxtLen
			idx := rdr.Idx
			start := rdr.Start

			// read next token
			for idx < txtlen {
				ch := text[idx]
				if ch == '<' {
					str := text[start:idx]
					idx++
					start = idx
					if IsNotJustWhitespace(str) {
						// update tokenizer
						rdr.Idx = idx
						rdr.Start = start
						// contents
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						if HasAmpOrNotASCII(str) {
							return html.UnescapeString(str), true
						} else {
							return str, true
						}
					}
				} else if ch == '>' {
					str := text[start:idx]
					idx++
					start = idx
					if IsNotJustWhitespace(str) {
						// update tokenizer
						rdr.Idx = idx
						rdr.Start = start
						// element
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						return str, false
					}
				} else {
					idx++
				}
			}

			// signal end of current block
			rdr.Text = ""

			// update tokenizer
			rdr.Idx = idx
			rdr.Start = start
		}

		return "", false
	}

	for {
		if rdr.EndToken != "" {
			// send end token for self-closing element
			name := rdr.EndToken
			rdr.EndToken = ""
			return IS_END_TAG, name, "", ""
		}

		// nextUnit removes < and > characters from element
		str, isContent := nextUnit()
		if str == "" {
			break
		}

		if isContent {
			return IS_CONTENT, "", "", str
		}

		if str[0] == '?' || str[0] == '!' {
			// skip ?xml, !DOCTYPE, !comment, and ?processing instructions
			continue
		}

		// look for leading or trailing / character
		isEnd := false
		isSelf := false

		if str[0] == '/' {
			isEnd = true
			str = str[1:]
		} else {
			strlen := len(str)
			if str[strlen-1] == '/' {
				isSelf = true
				str = str[0 : strlen-1]
			}
		}

		name := ""
		attributes := ""

		pos := strings.IndexByte(str, ' ')
		if pos > -1 {
			name = str[:pos]
			attributes = str[pos+1:]
		} else {
			name = str
		}

		if isEnd {
			return IS_END_TAG, name, "", ""
		}

		if isSelf {
			// queue end token
			rdr.EndToken = name
			// then proceed as a start token
		}

		return IS_START_TAG, name, attributes, ""
	}

	return IS_CLOSED, "", "", ""
}

// PARSE XML BLOCKS INTO STRINGS FROM <PATTERN> TO </PATTERN>

type Scanner struct {
	Pattern   string
	PatLength int
	CharSkip  [256]int
}

func NewScanner(pattern string) *Scanner {

	if pattern == "" {
		return nil
	}

	scr := &Scanner{Pattern: pattern}

	patlen := len(pattern)
	scr.PatLength = patlen

	// position of last character in pattern
	last := patlen - 1

	// initialize bad character displacement table
	for i := range scr.CharSkip {
		scr.CharSkip[i] = patlen
	}
	for i := 0; i < last; i++ {
		ch := pattern[i]
		scr.CharSkip[ch] = last - i
	}

	return scr
}

func IsAnElement(text string, lf, rt, mx int) bool {

	if (lf >= 0 && text[lf] == '<') || (lf > 0 && text[lf] == '/' && text[lf-1] == '<') {
		if (rt < mx && (text[rt] == '>' || text[rt] == ' ')) || (rt+1 < mx && text[rt] == '/' && text[rt+1] == '>') {
			return true
		}
	}

	return false
}

// modified Boyer-Moore-Horspool search function
func FindNextMatch(scr *Scanner, text string, offset int) (int, int, int) {

	if scr == nil || text == "" {
		return -1, -1, -1
	}

	// copy values into local variables for speed
	txtlen := len(text)
	pattern := scr.Pattern[:]
	patlen := scr.PatLength
	max := txtlen - patlen
	last := patlen - 1
	skip := scr.CharSkip[:]

	i := offset

	for i <= max {
		j := last
		k := i + last
		for j >= 0 && text[k] == pattern[j] {
			j--
			k--
		}
		// require match candidate to be element name, i.e., <pattern ... >, </pattern ... >, or <pattern ... />
		if j < 0 && IsAnElement(text, i-1, i+patlen, txtlen) {
			// find positions of flanking brackets
			lf := i - 1
			for lf > 0 && text[lf] != '<' {
				lf--
			}
			rt := i + patlen
			for rt < txtlen && text[rt] != '>' {
				rt++
			}
			return i + 1, lf, rt + 1
		}
		// find character in text above last character in pattern
		ch := text[i+last]
		// displacement table can shift pattern by one or more positions
		i += skip[ch]
	}

	return -1, -1, -1
}

// function to find next element with pattern name
func NextPattern(scr *Scanner, text string, pos int) (int, int, int) {

	if scr == nil || text == "" {
		return IS_NOT_A_TAG, 0, 0
	}

	prev := pos

	for {
		next, start, stop := FindNextMatch(scr, text, prev)
		if next < 0 {
			return IS_NOT_A_TAG, 0, 0
		}

		prev = next + 1

		if text[start+1] == '/' {
			return IS_END_TAG, stop, prev
		} else if text[stop-2] == '/' {
			return IS_SELF_TAG, start, prev
		} else {
			return IS_START_TAG, start, prev
		}
	}
}

// master partitioning function
func PartitionPattern(pat, star string, rdr *XmlReader, proc func(int, string)) {

	if pat == "" || rdr == nil || proc == nil {
		return
	}

	// handle -pattern Object construct
	doNormal := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator bytes.Buffer

		match := 0
		pos := 0
		next := 0

		rec := 0

		scr := NewScanner(pat)
		if scr == nil {
			return
		}

		for {

			begin = 0
			next = 0

			line = NextBlock(rdr)
			if line == "" {
				return
			}

			for {
				match, pos, next = NextPattern(scr, line, next)
				if match == IS_START_TAG {
					if level == 0 {
						inPattern = true
						begin = pos
					}
					level++
				} else if match == IS_END_TAG {
					level--
					if level == 0 {
						inPattern = false
						accumulator.WriteString(line[begin:pos])
						// read and process one -pattern object at a time
						str := accumulator.String()
						if str != "" {
							rec++
							proc(rec, str)
						}
						// reset accumulator
						accumulator.Reset()
					}
				} else if match == IS_SELF_TAG {
				} else {
					if inPattern {
						accumulator.WriteString(line[begin:])
					}
					break
				}
			}
		}
	}

	// handle -pattern Parent/* construct
	doStar := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator bytes.Buffer

		match := 0
		pos := 0
		next := 0

		rec := 0

		scr := NewScanner(pat)
		if scr == nil {
			return
		}

		last := pat

		// read to first <pattern> element
		for {

			next = 0

			line = NextBlock(rdr)
			if line == "" {
				break
			}

			match, pos, next = NextPattern(scr, line, next)
			if match == IS_START_TAG {
				break
			}
		}

		if match != IS_START_TAG {
			return
		}

		// function to find next element in XML
		nextElement := func(text string, pos int) string {

			txtlen := len(text)

			tag := ""
			for i := pos; i < txtlen; i++ {
				if text[i] == '<' {
					tag = text[i+1:]
					break
				}
			}
			if tag == "" {
				return ""
			}
			if tag[0] == '/' {
				return ""
			}
			for i, ch := range tag {
				if ch == '>' || ch == ' ' || ch == '/' {
					return tag[0:i]
				}
			}

			return ""
		}

		// read and process heterogeneous objects immediately below <pattern> parent
		for {
			tag := nextElement(line, next)
			if tag == "" {

				begin = 0
				next = 0

				line = NextBlock(rdr)
				if line == "" {
					break
				}

				tag = nextElement(line, next)
			}
			if tag == "" {
				return
			}
			if tag[0] == '/' {
				return
			}
			if tag != last {
				scr = NewScanner(tag)
				if scr == nil {
					return
				}
				last = tag
			}

			for {
				match, pos, next = NextPattern(scr, line, next)
				if match == IS_START_TAG {
					if level == 0 {
						inPattern = true
						begin = pos
					}
					level++
				} else if match == IS_END_TAG {
					level--
					if level == 0 {
						inPattern = false
						accumulator.WriteString(line[begin:pos])
						// read and process one -pattern/* object at a time
						str := accumulator.String()
						if str != "" {
							rec++
							proc(rec, str)
						}
						// reset accumulator
						accumulator.Reset()
						break
					}
				} else {
					if inPattern {
						accumulator.WriteString(line[begin:])
					}

					begin = 0
					next = 0

					line = NextBlock(rdr)
					if line == "" {
						break
					}
				}
			}
		}
	}

	// call appropriate handler
	if star == "" {
		doNormal()
	} else if star == "*" {
		doStar()
	}
}

// XML CLEANUP AND FORMATTING FUNCTION

// reformat XML for ease of reading
func ProcessFormat(doCompress, doCleanup bool) {

	in := bufio.NewReader(os.Stdin)
	blk := NewXmlReader(in, doCompress, doCleanup)
	rdr := NewTokenizer(blk, "")
	if rdr == nil {
		return
	}

	// function returns true for contents, false for element, also filters out <? and <! sections
	nextUnit := func() (string, bool) {

		skipTo := ""

		for {
			if rdr.Text == "" && rdr.Reader != nil {
				// get next block
				rdr.Text = NextBlock(rdr.Reader)
				rdr.TxtLen = len(rdr.Text)
				rdr.Idx = 0
				rdr.Start = 0
			}

			if rdr.Text == "" {
				// if closed, break out of loop
				break
			}

			// copy current tokenizer values into local variables for speed
			text := rdr.Text[:]
			txtlen := rdr.TxtLen
			idx := rdr.Idx
			start := rdr.Start

			// read next token
			for idx < txtlen {

				if skipTo != "" {
					found := strings.Index(text[idx:], skipTo)
					if found > -1 {
						found += idx
						// skip past end of region
						idx = found + len(skipTo)
						start = idx
						skipTo = ""
					} else {
						// continue looking on next block
						break
					}
				}

				ch := text[idx]
				if ch == '<' {
					str := text[start:idx]
					start = idx
					idx++
					if IsNotJustWhitespace(str) {
						// check for start of region to remove
						if strings.HasPrefix(str, "<!") {
							if strings.HasPrefix(str, "<!--") {
								skipTo = "-->"
							} else if strings.HasPrefix(str, "<![CDATA[") {
								skipTo = "]]>"
							}
							continue
						} else if strings.HasPrefix(str, "<?") {
							continue
						}
						// update tokenizer
						rdr.Idx = idx
						rdr.Start = start
						// contents
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						if HasAmpOrNotASCII(str) {
							return html.UnescapeString(str), true
						} else {
							return str, true
						}
					}
				} else if ch == '>' {
					idx++
					str := text[start:idx]
					start = idx
					if IsNotJustWhitespace(str) {
						// check for start of region to remove
						if strings.HasPrefix(str, "<!") {
							if strings.HasPrefix(str, "<!--") {
								skipTo = "-->"
							} else if strings.HasPrefix(str, "<![CDATA[") {
								skipTo = "]]>"
							}
							continue
						} else if strings.HasPrefix(str, "<?") {
							continue
						}
						// update tokenizer
						rdr.Idx = idx
						rdr.Start = start
						// remove < and >
						str = str[1 : len(str)-1]
						// element
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						return str, false
					}
				} else {
					idx++
				}
			}

			// signal end of current block
			rdr.Text = ""

			// update tokenizer
			rdr.Idx = idx
			rdr.Start = start
		}

		return "", false
	}

	// return class, element name, attributes, content text
	nextFormatToken := func() (int, string, string, string) {

		for {
			if rdr.EndToken != "" {
				// send end token for self-closing element
				name := rdr.EndToken
				rdr.EndToken = ""
				return IS_END_TAG, name, "", ""
			}

			// nextUnit removes < and > characters from element
			str, isContent := nextUnit()
			if str == "" {
				break
			}

			if isContent {
				return IS_CONTENT, "", "", str
			}

			if str[0] == '?' || str[0] == '!' {
				// <? and <! constructs should already be filtered out by -format version of nextUnit
				continue
			}

			// look for leading or trailing / character
			isEnd := false
			isSelf := false

			if str[0] == '/' {
				isEnd = true
				str = str[1:]
			} else {
				strlen := len(str)
				if str[strlen-1] == '/' {
					isSelf = true
					str = str[0 : strlen-1]
				}
			}

			name := ""
			attributes := ""

			pos := strings.IndexByte(str, ' ')
			if pos > -1 {
				name = str[:pos]
				attributes = str[pos+1:]
			} else {
				name = str
			}

			if isEnd {
				return IS_END_TAG, name, "", ""
			}

			if isSelf {
				// queue end token
				rdr.EndToken = name
				// then proceed as a start token
			}

			return IS_START_TAG, name, attributes, ""
		}

		return IS_CLOSED, "", "", ""
	}

	const (
		NOTSET = iota
		START
		END
		CHAR
		OTHER
	)

	// array to speed up indentation
	indentSpaces := []string{
		"",
		"  ",
		"    ",
		"      ",
		"        ",
		"          ",
		"            ",
		"              ",
		"                ",
		"                  ",
	}

	indent := 0

	// parent used to detect first start tag, will place in doctype line
	parent := ""

	var buffer bytes.Buffer
	count := 0

	status := NOTSET

	// delay printing right bracket of start tag to support self-closing tag style
	needsRightBracket := ""

	// delay printing start tag if no attributes, suppress empty start-end pair if followed by end
	justStartName := ""
	justStartIndent := 0

	// function to indent a specified number of spaces
	doIndent := func(indt int) {
		i := indt
		for i > 9 {
			buffer.WriteString("                    ")
			i -= 10
		}
		if i < 0 {
			return
		}
		buffer.WriteString(indentSpaces[i])
	}

	// function to handle delayed start tag
	doDelayedName := func() {
		if needsRightBracket != "" {
			buffer.WriteString(">")
			needsRightBracket = ""
		}
		if justStartName != "" {
			doIndent(justStartIndent)
			buffer.WriteString("<")
			buffer.WriteString(justStartName)
			buffer.WriteString(">")
			justStartName = ""
		}
	}

	for {
		cat, name, attributes, contents := nextFormatToken()

		switch cat {
		case IS_START_TAG:
			doDelayedName()
			if status == START {
				buffer.WriteString("\n")
			}
			// remove internal copies of </parent><parent> tags
			if parent != "" && name == parent && indent == 1 {
				continue
			}

			// detect first start tag, print xml and doctype parent
			if indent == 0 && parent == "" {
				parent = name
				fmt.Printf("<?xml version=\"1.0\"?>\n")
				fmt.Printf("<!DOCTYPE %s>\n", parent)
				// now filtering internal </parent><parent> tags, so queue printing of closing tag
				defer fmt.Printf("</%s>\n", parent)
				// already past </parent><parent> test, so opening tag will print normally
			}

			if IsNotJustWhitespace(attributes) {
				doIndent(indent)

				buffer.WriteString("<")
				buffer.WriteString(name)

				attr := strings.TrimSpace(attributes)
				attr = CompressRunsOfSpaces(attr)
				buffer.WriteString(" ")
				buffer.WriteString(attr)
				needsRightBracket = name
			} else {
				justStartName = name
				justStartIndent = indent
			}

			indent++

			status = START
		case IS_END_TAG:
			// if end immediately follows start, turn into self-closing tag if there were attributes, otherwise suppress empty tag
			if needsRightBracket != "" {
				if status == START && name == needsRightBracket {
					// end immediately follows start, produce self-closing tag
					buffer.WriteString("/>\n")
					needsRightBracket = ""
					indent--
					status = END
					break
				}
				buffer.WriteString(">")
				needsRightBracket = ""
			}
			if justStartName != "" {
				if status == START && name == justStartName {
					// end immediately follows delayed start with no attributes, suppress
					justStartName = ""
					indent--
					status = END
					break
				}
				doIndent(justStartIndent)
				buffer.WriteString("<")
				buffer.WriteString(justStartName)
				buffer.WriteString(">")
				justStartName = ""
			}

			// remove internal copies of </parent><parent> tags
			if parent != "" && name == parent && indent == 1 {
				continue
			}
			indent--
			if status == CHAR {
				buffer.WriteString("</")
				buffer.WriteString(name)
				buffer.WriteString(">\n")
			} else if status == START {
				buffer.WriteString("</")
				buffer.WriteString(name)
				buffer.WriteString(">\n")
			} else {
				doIndent(indent)

				buffer.WriteString("</")
				buffer.WriteString(name)
				buffer.WriteString(">\n")
			}
			status = END
		case IS_CONTENT:
			doDelayedName()
			if len(contents) > 0 && IsNotJustWhitespace(contents) {
				contents = html.EscapeString(contents)
				buffer.WriteString(contents)
				status = CHAR
			}
		case IS_CLOSED:
			doDelayedName()
			txt := buffer.String()
			if txt != "" {
				// print final buffer
				fmt.Printf("%s", txt)
			}
			return
		default:
			doDelayedName()
			status = OTHER
		}

		count++
		if count > 1000 {
			count = 0
			txt := buffer.String()
			if txt != "" {
				// print current buffered output
				fmt.Printf("%s", txt)
			}
			buffer.Reset()
		}
	}
}

// MISCELLANEOUS XML PRESENTATION AND VALIDATION FUNCTIONS

// display outline of XML structure
func ProcessOutline(doCompress, doCleanup bool) {

	in := bufio.NewReader(os.Stdin)
	blk := NewXmlReader(in, doCompress, doCleanup)
	rdr := NewTokenizer(blk, "")
	if rdr == nil {
		return
	}

	indent := 0

	for {
		cat, name, _, _ := NextToken(rdr)
		switch cat {
		case IS_START_TAG:
			if name == "eSummaryResult" ||
				name == "eLinkResult" ||
				name == "eInfoResult" ||
				name == "PubmedArticleSet" ||
				name == "DocumentSummarySet" ||
				name == "INSDSet" ||
				name == "Entrezgene-Set" ||
				name == "TaxaSet" {
				continue
			}
			for i := 0; i < indent; i++ {
				fmt.Printf("  ")
			}
			fmt.Printf("%s\n", name)
			indent++
		case IS_END_TAG:
			indent--
		case IS_CLOSED:
			return
		default:
		}
	}
}

// display paths to XML elements
func SynopsisLevel(parent string, rdr *Tokenizer) bool {

	for {
		cat, name, _, _ := NextToken(rdr)
		switch cat {
		case IS_START_TAG:
			if name == "eSummaryResult" ||
				name == "eLinkResult" ||
				name == "eInfoResult" ||
				name == "PubmedArticleSet" ||
				name == "DocumentSummarySet" ||
				name == "INSDSet" ||
				name == "Entrezgene-Set" ||
				name == "TaxaSet" {
				continue
			}
			if parent != "" {
				fmt.Printf("%s/", parent)
			}
			fmt.Printf("%s\n", name)
			path := parent
			if path != "" {
				path += "/"
			}
			path += name
			is_closed := SynopsisLevel(path, rdr)
			if is_closed {
				return true
			}
		case IS_END_TAG:
			// break recursion
			return false
		case IS_CLOSED:
			return true
		default:
		}
	}
}

func ProcessSynopsis(doCompress, doCleanup bool) {

	in := bufio.NewReader(os.Stdin)
	blk := NewXmlReader(in, doCompress, doCleanup)
	rdr := NewTokenizer(blk, "")
	if rdr == nil {
		return
	}

	for {
		is_closed := SynopsisLevel("", rdr)
		if is_closed {
			return
		}
	}
}

// verify integrity of XML object nesting
func VerifyLevel(parent string, level int, rdr *Tokenizer) {

	// function returns true for contents, false for element, counts newline characters
	nextUnit := func() (string, bool) {

		for {
			if rdr.Text == "" && rdr.Reader != nil {
				// get next block
				rdr.Text = NextBlock(rdr.Reader)
				rdr.TxtLen = len(rdr.Text)
				rdr.Idx = 0
				rdr.Start = 0
			}

			if rdr.Text == "" {
				// if closed, break out of loop
				break
			}

			// copy current tokenizer values into local variables for speed
			text := rdr.Text[:]
			txtlen := rdr.TxtLen
			idx := rdr.Idx
			start := rdr.Start

			// read next token
			for idx < txtlen {
				ch := text[idx]
				if ch == '<' {
					str := text[start:idx]
					idx++
					start = idx
					if IsNotJustWhitespace(str) {
						// update tokenizer
						rdr.Idx = idx
						rdr.Start = start
						// contents
						// no need to return actual contents for XML element name balance test
						return "", true
					}
				} else if ch == '>' {
					str := text[start:idx]
					idx++
					start = idx
					if IsNotJustWhitespace(str) {
						// update tokenizer
						rdr.Idx = idx
						rdr.Start = start
						// element
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						return str, false
					}
				} else if ch == '\n' {
					// expect Unix style single newline character
					rdr.Line++
					idx++
				} else {
					idx++
				}
			}

			// signal end of current block
			rdr.Text = ""

			// update tokenizer
			rdr.Idx = idx
			rdr.Start = start
		}

		return "", false
	}

	// for XML integrity check, return class and element name, keeping track of line number
	nextVerifyToken := func() (int, string) {

		for {
			if rdr.EndToken != "" {
				// send end token for self-closing element
				name := rdr.EndToken
				rdr.EndToken = ""
				return IS_END_TAG, name
			}

			// nextUnit removes < and > characters from element
			str, isContent := nextUnit()

			if isContent {
				return IS_CONTENT, ""
			}

			// if not content, need element name
			if str == "" {
				break
			}

			if str[0] == '?' || str[0] == '!' {
				// skip ?xml, !DOCTYPE, !comment, and ?processing instructions
				continue
			}

			// look for leading or trailing / character
			isEnd := false
			isSelf := false

			if str[0] == '/' {
				isEnd = true
				str = str[1:]
			} else {
				strlen := len(str)
				if str[strlen-1] == '/' {
					isSelf = true
					str = str[0 : strlen-1]
				}
			}

			// isolate element name, not checking attribute validity
			pos := strings.IndexByte(str, ' ')
			if pos > -1 {
				str = str[:pos]
			}

			if isEnd {
				return IS_END_TAG, str
			}

			if isSelf {
				// queue end token
				rdr.EndToken = str
				// then proceed as a start token
			}

			return IS_START_TAG, str
		}

		return IS_CLOSED, ""
	}

	const (
		NOTSET = iota
		START
		END
		CHAR
		OTHER
	)

	status := START
	for {
		cat, name := nextVerifyToken()
		line := rdr.Line

		switch cat {
		case IS_START_TAG:
			if status == CHAR {
				fmt.Printf("<%s> not expected after contents, line %d\n", name, line)
			}
			VerifyLevel(name, level+1, rdr)
			// returns here after recursion
			status = END
		case IS_END_TAG:
			if parent != name && parent != "" {
				fmt.Printf("Expected </%s>, found </%s>, line %d\n", parent, name, line)
			}
			if level < 1 {
				fmt.Printf("Unexpected </%s> at end of XML, line %d\n", name, line)
			}
			// break recursion
			return
		case IS_CONTENT:
			if status != START {
				fmt.Printf("Contents not expected before </%s>, line %d\n", parent, line)
			}
			status = CHAR
		case IS_CLOSED:
			if level > 0 {
				fmt.Printf("Unexpected end of data, line %d\n", line)
			}
			return
		default:
			status = OTHER
		}
	}
}

func ProcessVerify(doCompress, doCleanup bool) {

	in := bufio.NewReader(os.Stdin)
	blk := NewXmlReader(in, doCompress, doCleanup)
	rdr := NewTokenizer(blk, "")
	if rdr == nil {
		return
	}

	VerifyLevel("", 0, rdr)
}

// INSDSEQ EXTRACTION COMMAND GENERATOR

// e.g., xtract -insd complete mat_peptide "%peptide" product peptide

func CheckAgainstVocabulary(str, objtype string, arry []string) {

	if str == "" || arry == nil {
		return
	}

	// skip past pound, percent, or caret character at beginning of string
	if len(str) > 1 {
		switch str[0] {
		case '#', '%', '^':
			str = str[1:]
		default:
		}
	}

	for _, txt := range arry {
		if str == txt {
			return
		}
		if strings.ToUpper(str) == strings.ToUpper(txt) {
			fmt.Fprintf(os.Stderr, "\nIncorrect capitalization of '%s' %s, change to '%s'\n", str, objtype, txt)
			os.Exit(1)
		}
	}

	fmt.Fprintf(os.Stderr, "\nItem '%s' is not a legal -insd %s\n", str, objtype)
	os.Exit(1)
}

func ProcessINSD(args []string, isPipe bool) []string {

	// legal GenBank / GenPept / RefSeq features

	features := []string{
		"-10_signal",
		"-35_signal",
		"3'clip",
		"3'UTR",
		"5'clip",
		"5'UTR",
		"allele",
		"assembly_gap",
		"attenuator",
		"C_region",
		"CAAT_signal",
		"CDS",
		"centromere",
		"conflict",
		"D_segment",
		"D-loop",
		"enhancer",
		"exon",
		"gap",
		"GC_signal",
		"gene",
		"iDNA",
		"intron",
		"J_segment",
		"LTR",
		"mat_peptide",
		"misc_binding",
		"misc_difference",
		"misc_feature",
		"misc_recomb",
		"misc_RNA",
		"misc_signal",
		"misc_structure",
		"mobile_element",
		"modified_base",
		"mRNA",
		"mutation",
		"N_region",
		"ncRNA",
		"old_sequence",
		"operon",
		"oriT",
		"polyA_signal",
		"polyA_site",
		"precursor_RNA",
		"prim_transcript",
		"primer_bind",
		"promoter",
		"protein_bind",
		"RBS",
		"regulatory",
		"rep_origin",
		"repeat_region",
		"repeat_unit",
		"rRNA",
		"S_region",
		"satellite",
		"scRNA",
		"sig_peptide",
		"snoRNA",
		"snRNA",
		"source",
		"stem_loop",
		"STS",
		"TATA_signal",
		"telomere",
		"terminator",
		"tmRNA",
		"transit_peptide",
		"tRNA",
		"unsure",
		"V_region",
		"V_segment",
		"variation",
	}

	// legal GenBank / GenPept / RefSeq qualifiers

	qualifiers := []string{
		"allele",
		"altitude",
		"anticodon",
		"artificial_location",
		"bio_material",
		"bound_moiety",
		"calculated_mol_wt",
		"cell_line",
		"cell_type",
		"chloroplast",
		"chromoplast",
		"chromosome",
		"citation",
		"clone_lib",
		"clone",
		"codon_start",
		"codon",
		"collected_by",
		"collection_date",
		"compare",
		"cons_splice",
		"country",
		"cultivar",
		"culture_collection",
		"cyanelle",
		"db_xref",
		"dev_stage",
		"direction",
		"EC_number",
		"ecotype",
		"environmental_sample",
		"estimated_length",
		"evidence",
		"exception",
		"experiment",
		"focus",
		"frequency",
		"function",
		"gap_type",
		"gdb_xref",
		"gene_synonym",
		"gene",
		"germline",
		"haplogroup",
		"haplotype",
		"identified_by",
		"inference",
		"insertion_seq",
		"isolate",
		"isolation_source",
		"kinetoplast",
		"lab_host",
		"label",
		"lat_lon",
		"linkage_evidence",
		"locus_tag",
		"macronuclear",
		"map",
		"mating_type",
		"metagenome_source",
		"metagenomic",
		"mitochondrion",
		"mobile_element_type",
		"mobile_element",
		"mod_base",
		"mol_type",
		"ncRNA_class",
		"non_functional",
		"note",
		"number",
		"old_locus_tag",
		"operon",
		"organelle",
		"organism",
		"partial",
		"PCR_conditions",
		"PCR_primers",
		"peptide",
		"phenotype",
		"plasmid",
		"pop_variant",
		"product",
		"protein_id",
		"proviral",
		"pseudo",
		"pseudogene",
		"rearranged",
		"regulatory_class",
		"replace",
		"ribosomal_slippage",
		"rpt_family",
		"rpt_type",
		"rpt_unit_range",
		"rpt_unit_seq",
		"rpt_unit",
		"satellite",
		"segment",
		"sequenced_mol",
		"serotype",
		"serovar",
		"sex",
		"specific_host",
		"specimen_voucher",
		"standard_name",
		"strain",
		"sub_clone",
		"sub_species",
		"sub_strain",
		"tag_peptide",
		"tissue_lib",
		"tissue_type",
		"trans_splicing",
		"transcript_id",
		"transcription",
		"transgenic",
		"transl_except",
		"transl_table",
		"translation",
		"transposon",
		"type_material",
		"UniProtKB_evidence",
		"usedin",
		"variety",
		"virion",
	}

	// legal INSDSeq XML fields

	insdtags := []string{
		"INSDAltSeqData_items",
		"INSDAltSeqData",
		"INSDAltSeqItem_first-accn",
		"INSDAltSeqItem_gap-comment",
		"INSDAltSeqItem_gap-length",
		"INSDAltSeqItem_gap-linkage",
		"INSDAltSeqItem_gap-type",
		"INSDAltSeqItem_interval",
		"INSDAltSeqItem_isgap",
		"INSDAltSeqItem_last-accn",
		"INSDAltSeqItem_value",
		"INSDAltSeqItem",
		"INSDAuthor",
		"INSDComment_paragraphs",
		"INSDComment_type",
		"INSDComment",
		"INSDCommentParagraph",
		"INSDFeature_intervals",
		"INSDFeature_key",
		"INSDFeature_location",
		"INSDFeature_operator",
		"INSDFeature_partial3",
		"INSDFeature_partial5",
		"INSDFeature_quals",
		"INSDFeature_xrefs",
		"INSDFeature",
		"INSDFeatureSet_annot-source",
		"INSDFeatureSet_features",
		"INSDFeatureSet",
		"INSDInterval_accession",
		"INSDInterval_from",
		"INSDInterval_interbp",
		"INSDInterval_iscomp",
		"INSDInterval_point",
		"INSDInterval_to",
		"INSDInterval",
		"INSDKeyword",
		"INSDQualifier_name",
		"INSDQualifier_value",
		"INSDQualifier",
		"INSDReference_authors",
		"INSDReference_consortium",
		"INSDReference_journal",
		"INSDReference_position",
		"INSDReference_pubmed",
		"INSDReference_reference",
		"INSDReference_remark",
		"INSDReference_title",
		"INSDReference_xref",
		"INSDReference",
		"INSDSecondary-accn",
		"INSDSeq_accession-version",
		"INSDSeq_alt-seq",
		"INSDSeq_comment-set",
		"INSDSeq_comment",
		"INSDSeq_contig",
		"INSDSeq_create-date",
		"INSDSeq_create-release",
		"INSDSeq_database-reference",
		"INSDSeq_definition",
		"INSDSeq_division",
		"INSDSeq_entry-version",
		"INSDSeq_feature-set",
		"INSDSeq_feature-table",
		"INSDSeq_keywords",
		"INSDSeq_length",
		"INSDSeq_locus",
		"INSDSeq_moltype",
		"INSDSeq_organism",
		"INSDSeq_other-seqids",
		"INSDSeq_primary-accession",
		"INSDSeq_primary",
		"INSDSeq_project",
		"INSDSeq_references",
		"INSDSeq_secondary-accessions",
		"INSDSeq_segment",
		"INSDSeq_sequence",
		"INSDSeq_source-db",
		"INSDSeq_source",
		"INSDSeq_strandedness",
		"INSDSeq_struc-comments",
		"INSDSeq_taxonomy",
		"INSDSeq_topology",
		"INSDSeq_update-date",
		"INSDSeq_update-release",
		"INSDSeq_xrefs",
		"INSDSeq",
		"INSDSeqid",
		"INSDSet",
		"INSDStrucComment_items",
		"INSDStrucComment_name",
		"INSDStrucComment",
		"INSDStrucCommentItem_tag",
		"INSDStrucCommentItem_url",
		"INSDStrucCommentItem_value",
		"INSDStrucCommentItem",
		"INSDXref_dbname",
		"INSDXref_id",
		"INSDXref",
	}

	var acc []string

	max := len(args)
	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nInsufficient command-line arguments supplied to xtract -insd\n")
		os.Exit(1)
	}

	acc = append(acc, "-pattern", "INSDSeq", "-ACCN", "INSDSeq_accession-version")
	printAccn := true

	// collect descriptors

	if strings.HasPrefix(args[0], "INSD") {

		if isPipe {
			acc = append(acc, "-clr", "-pfx", "\\n", "-element", "&ACCN")
			acc = append(acc, "-group", "INSDSeq", "-sep", "|", "-element")
		} else {
			acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-element", "\"&ACCN\"")
			acc = append(acc, "-group", "INSDSeq", "-sep", "\"|\"", "-element")
		}
		printAccn = false

		for {
			if len(args) < 1 {
				return acc
			}
			str := args[0]
			if !strings.HasPrefix(args[0], "INSD") {
				break
			}
			CheckAgainstVocabulary(str, "element", insdtags)
			acc = append(acc, str)
			args = args[1:]
		}

	} else if strings.HasPrefix(strings.ToUpper(args[0]), "INSD") {

		// report capitalization or vocabulary failure
		CheckAgainstVocabulary(args[0], "element", insdtags)

		// program should not get to this point, but warn and exit anyway
		fmt.Fprintf(os.Stderr, "\nItem '%s' is not a legal -insd %s\n", args[0], "element")
		os.Exit(1)
	}

	// collect qualifiers

	partial := false
	complete := false

	if args[0] == "+" || args[0] == "complete" {
		complete = true
		args = args[1:]
		max--
	} else if args[0] == "-" || args[0] == "partial" {
		partial = true
		args = args[1:]
		max--
	}

	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nNo feature key supplied to xtract -insd\n")
		os.Exit(1)
	}

	acc = append(acc, "-group", "INSDFeature")

	// limit to designated features

	feature := args[0]

	fcmd := "-match"

	// can specify multiple features separated by plus sign (e.g., CDS+mRNA) or comma (e.g., CDS,mRNA)
	plus := strings.Split(feature, "+")
	for _, pls := range plus {
		comma := strings.Split(pls, ",")
		for _, cma := range comma {

			CheckAgainstVocabulary(cma, "feature", features)
			acc = append(acc, fcmd)
			ft := fmt.Sprintf("INSDFeature_key:%s", cma)
			acc = append(acc, ft)

			fcmd = "-or"
		}
	}

	if max < 2 {
		// still need at least one qualifier even on legal feature
		fmt.Fprintf(os.Stderr, "\nFeature '%s' must be followed by at least one qualifier\n", feature)
		os.Exit(1)
	}

	args = args[1:]

	if complete {
		acc = append(acc, "-avoid", "INSDFeature_partial5", "-and", "INSDFeature_partial3")
	} else if partial {
		acc = append(acc, "-match", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
	}

	if printAccn {
		if isPipe {
			acc = append(acc, "-clr", "-pfx", "\\n", "-element", "&ACCN")
		} else {
			acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-element", "\"&ACCN\"")
		}
	}

	for _, str := range args {
		if strings.HasPrefix(str, "INSD") {

			CheckAgainstVocabulary(str, "element", insdtags)
			if isPipe {
				acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
			} else {
				acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
			}
			acc = append(acc, str)

		} else if strings.HasPrefix(str, "#INSD") {

			CheckAgainstVocabulary(str, "element", insdtags)
			if isPipe {
				acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
				acc = append(acc, str)
			} else {
				acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
				ql := fmt.Sprintf("\"%s\"", str)
				acc = append(acc, ql)
			}

		} else if strings.HasPrefix(strings.ToUpper(str), "#INSD") || strings.HasPrefix(strings.ToUpper(str), "#INSD") {

			// report capitalization or vocabulary failure
			CheckAgainstVocabulary(str, "element", insdtags)

		} else {

			acc = append(acc, "-block", "INSDQualifier")

			CheckAgainstVocabulary(str, "qualifier", qualifiers)
			if len(str) > 2 && str[0] == '%' {
				acc = append(acc, "-match")
				ql := fmt.Sprintf("INSDQualifier_name:%s", str[1:])
				acc = append(acc, ql)
				if isPipe {
					acc = append(acc, "-element", "%INSDQualifier_value")
				} else {
					acc = append(acc, "-element", "\"%INSDQualifier_value\"")
				}
			} else {
				acc = append(acc, "-match")
				ql := fmt.Sprintf("INSDQualifier_name:%s", str)
				acc = append(acc, ql)
				acc = append(acc, "-element", "INSDQualifier_value")
			}
		}
	}

	return acc
}

// HYDRA CITATION MATCHER COMMAND GENERATOR

// acceptable scores are 0.8 or higher, exact match on "1" rejects low value in scientific notation with minus sign present
func ProcessHydra(isPipe bool) []string {

	var acc []string

	acc = append(acc, "-pattern", "Id")
	acc = append(acc, "-match", "@score", "-equals", "1")
	acc = append(acc, "-or", "@score", "-starts-with", "0.9")
	acc = append(acc, "-or", "@score", "-starts-with", "0.8")
	acc = append(acc, "-element", "Id")

	return acc
}

// PARSE COMMAND-LINE ARGUMENTS

/*
xtract -pattern Field -match "IsDate:Y" -and "IsHidden:N" -pfx "[" -sep "]\t" -element Name,FullName

<Block>
  <Visit> Field </Visit>
  <Conditions> -match IsDate:Y -and IsHidden:N <Conditions>
  <Subtasks>
    <Block>
      <Commands> -pfx [ -sep ]\t -element Name,FullName <Commands>
    </Block>
  </Subtasks>
</Block>

Block.Conditions and Block.Commands now parsed from Block.Tests and Block.Actions, respectively
*/

const (
	_ = iota
	UNIT
	SUBSET
	SECTION
	BLOCK
	BRANCH
	GROUP
	DIVISION
	PATTERN
)

// different names of exploration control arguments allow multiple levels of nested "for" loops in linear command line
// (capitalized versions for backward-compatibility with original Perl implementation handling of recursive definitions)
var (
	lcname = []string{"", "-unit", "-subset", "-section", "-block", "-branch", "-group", "-division", "-pattern"}

	ucname = []string{"", "-Unit", "-Subset", "-Section", "-Block", "-Branch", "-Group", "-Division", "-Pattern"}
)

func ParseCommands(parent *Block, startLevel int) {

	// function to find next highest level exploration argument
	findNextLevel := func(args []string, level int) (int, string, string) {

		if len(args) > 1 {

			for level > 0 {

				lctag := lcname[level]
				uctag := ucname[level]

				for _, txt := range args {
					if txt == lctag || txt == uctag {
						return level, lctag, uctag
					}
				}

				level--
			}
		}

		return 0, "", ""
	}

	arguments := parent.Actions

	level, lctag, uctag := findNextLevel(arguments, startLevel)

	if level < 1 {

		if parent.Visit != "" && len(arguments) > 0 {

			blk := &Block{Actions: arguments}
			parent.Subtasks = append(parent.Subtasks, blk)
			parent.Actions = nil
		}

		// break recursion
		return
	}

	// function to group arguments at a given exploration level
	subsetCommands := func(args []string, tagMatch bool) *Block {

		if !tagMatch {
			return &Block{Actions: args}
		}

		max := len(args)

		visit := ""

		// extract name of object to visit
		if max > 1 {
			visit = args[1]
			args = args[2:]
			max -= 2
		}

		if max < 1 {
			return &Block{Visit: visit}
		}

		position := ""

		// extract position argument
		if max > 1 && args[0] == "-position" {
			position = args[1]
			args = args[2:]
			max -= 2
		}

		if max < 1 {
			return &Block{Visit: visit, Position: position}
		}

		partition := 0
		for cur, str := range args {

			// record point between conditionals and remaining commands
			partition = cur + 1

			// skip if not a command
			if len(str) < 1 || str[0] != '-' {
				continue
			}

			if str == "-match" ||
				str == "-avoid" ||
				str == "-and" ||
				str == "-or" ||
				str == "-equals" ||
				str == "-contains" ||
				str == "-starts-with" ||
				str == "-ends-with" ||
				str == "-gt" ||
				str == "-ge" ||
				str == "-lt" ||
				str == "-le" ||
				str == "-eq" ||
				str == "-ne" {
				continue
			}

			// -element or other execution argument
			partition = cur

			// break out of loop
			break
		}

		// separate conditional and execution arguments
		conditions := args[0:partition]
		commands := args[partition:]

		return &Block{Visit: visit, Position: position, Tests: conditions, Actions: commands}
	}

	cur := 0

	tagMatch := false

	// search for positions of current exploration command

	for idx, txt := range arguments {
		if txt == lctag || txt == uctag {
			if idx == 0 {
				tagMatch = true
				continue
			}

			blk := subsetCommands(arguments[cur:idx], tagMatch)
			ParseCommands(blk, level-1)
			parent.Subtasks = append(parent.Subtasks, blk)

			cur = idx
			tagMatch = true
		}
	}

	if cur < len(arguments) {
		blk := subsetCommands(arguments[cur:], tagMatch)
		ParseCommands(blk, level-1)
		parent.Subtasks = append(parent.Subtasks, blk)
	}

	// clear execution arguments from parent after subsetting
	parent.Actions = nil
}

func ValidateCommands(cmds *Block) {

	if cmds == nil {
		return
	}

	// reality check on conditional arguments
	max := len(cmds.Tests)
	if max > 0 {

		// use inner exploration argument to match on contents if first selecting by position
		if cmds.Position != "" {
			fmt.Fprintf(os.Stderr, "\nCannot combine -position with -match or -avoid commands\n")
			os.Exit(1)
		}

		// check for missing -match (or -avoid) condition
		str := cmds.Tests[0]
		if str != "-match" && str != "-avoid" {
			fmt.Fprintf(os.Stderr, "\nMissing -match command before '%s'\n", str)
			os.Exit(1)
		}

		// check for missing argument after last -match (or -avoid) condition
		str = cmds.Tests[max-1]
		if len(str) > 0 && str[0] == '-' {
			fmt.Fprintf(os.Stderr, "\nItem missing after %s command\n", str)
			os.Exit(1)
		}
	}

	// conditionals should alternate between command and object/value
	expectDash := true
	last := ""
	for _, str := range cmds.Tests {
		if expectDash {
			if len(str) < 1 || str[0] != '-' {
				fmt.Fprintf(os.Stderr, "\nUnexpected '%s' argument after '%s'\n", str, last)
				os.Exit(1)
			}
			expectDash = false
		} else {
			if len(str) > 0 && str[0] == '-' {
				fmt.Fprintf(os.Stderr, "\nUnexpected '%s' command after '%s'\n", str, last)
				os.Exit(1)
			}
			expectDash = true
		}
		last = str
	}

	// reality check on execution arguments
	max = len(cmds.Actions)
	if max > 0 {

		// check for missing -element (or -first, etc.) command
		str := cmds.Actions[0]
		if len(str) < 1 || str[0] != '-' {
			fmt.Fprintf(os.Stderr, "\nMissing -element command before '%s'\n", str)
			os.Exit(1)
		}

		// check for -match after -element
		for _, str := range cmds.Actions {

			// skip if not a command
			if len(str) < 1 || str[0] != '-' {
				continue
			}

			if str == "-match" ||
				str == "-avoid" ||
				str == "-and" ||
				str == "-or" ||
				str == "-equals" ||
				str == "-contains" ||
				str == "-starts-with" ||
				str == "-ends-with" ||
				str == "-gt" ||
				str == "-ge" ||
				str == "-lt" ||
				str == "-le" ||
				str == "-eq" ||
				str == "-ne" {
				fmt.Fprintf(os.Stderr, "\nMisplaced %s command\n", str)
				os.Exit(1)
			}
		}

		// check for missing argument after last -element (or -first, etc.) command
		str = cmds.Actions[max-1]
		if len(str) > 0 && str[0] == '-' {
			if str == "-rst" {
				fmt.Fprintf(os.Stderr, "\nUnexpected position for %s command\n", str)
				os.Exit(1)
			} else if str != "-clr" {
				fmt.Fprintf(os.Stderr, "\nItem missing after %s command\n", str)
				os.Exit(1)
			}
		}
	}

	for _, sub := range cmds.Subtasks {
		ValidateCommands(sub)
	}
}

func ParseActions(cmds *Block) {

	arguments := cmds.Actions[:]

	max := len(arguments)
	if max < 1 {
		return
	}

	comm := make([]*Operation, 0, max)

	status := UNSET

	// function to parse next argument
	nextStatus := func(str string) int {
		status = ParseFlag(str)

		switch status {
		case VARIABLE:
			op := &Operation{Type: status, Value: str[1:]}
			comm = append(comm, op)
			status = VALUE
		case CLR, RST:
			op := &Operation{Type: status, Value: ""}
			comm = append(comm, op)
			status = UNSET
		case UNSET:
			fmt.Fprintf(os.Stderr, "\nNo -element before '%s'\n", str)
			os.Exit(1)
		case UNRECOGNIZED:
			fmt.Fprintf(os.Stderr, "\nUnrecognized argument '%s'\n", str)
			os.Exit(1)
		default:
		}

		return status
	}

	idx := 0

	// parse command strings into operation structure
	for idx < max {
		str := arguments[idx]
		idx++

		switch status {
		case UNSET:
			status = nextStatus(str)
		case ELEMENT, FIRST, LAST, ENCODE, SUM, MIN, MAX:
			for !strings.HasPrefix(str, "-") {
				op := &Operation{Type: status, Value: str}
				comm = append(comm, op)
				if idx >= max {
					break
				}
				str = arguments[idx]
				idx++
			}
			status = UNSET
			if idx < max {
				status = nextStatus(str)
			}
		case TAB, RET, PFX, SFX, SEP, LBL, PFC:
			op := &Operation{Type: status, Value: ConvertSlash(str)}
			comm = append(comm, op)
			status = UNSET
		case VARIABLE:
			op := &Operation{Type: status, Value: str[1:]}
			comm = append(comm, op)
			status = VALUE
		case VALUE:
			op := &Operation{Type: status, Value: str}
			comm = append(comm, op)
			status = UNSET
		case UNRECOGNIZED:
			fmt.Fprintf(os.Stderr, "\nUnrecognized argument '%s'\n", str)
			os.Exit(1)
		default:
		}
	}

	cmds.Commands = comm
}

func ParseTests(cmds *Block) {

	arguments := cmds.Tests[:]

	max := len(arguments)
	if max < 1 {
		return
	}

	cond := make([]*Operation, 0, max)

	status := UNSET

	idx := 0

	// parse command strings into operation structure
	for idx < max {
		str := arguments[idx]
		idx++

		switch status {
		case UNSET:
			status = ParseFlag(str)
		case MATCH, AVOID, AND, OR, EQUALS, CONTAINS, STARTSWITH, ENDSWITH, GT, GE, LT, LE, EQ, NE:
			op := &Operation{Type: status, Value: str}
			cond = append(cond, op)
			status = UNSET
		case UNRECOGNIZED:
			fmt.Fprintf(os.Stderr, "\nUnrecognized argument '%s'\n", str)
			os.Exit(1)
		default:
		}
	}

	cmds.Conditions = cond
}

func ResolveOperations(cmds *Block) {

	ParseActions(cmds)
	ParseTests(cmds)

	for _, sub := range cmds.Subtasks {
		ResolveOperations(sub)
	}
}

func ParseArguments(args []string) *Block {

	// parse nested exploration instruction from command-line arguments
	head := &Block{}

	for _, txt := range args {
		head.Actions = append(head.Actions, txt)
	}

	ParseCommands(head, PATTERN)

	// check for no -element or multiple -pattern commands
	noElement := true
	numPatterns := 0
	for _, txt := range args {
		if txt == "-element" ||
			txt == "-first" ||
			txt == "-last" ||
			txt == "-encode" ||
			txt == "-sum" ||
			txt == "-min" ||
			txt == "-max" {
			noElement = false
		}
		if txt == "-pattern" || txt == "-Pattern" {
			numPatterns++
		}
	}

	if numPatterns > 1 {
		fmt.Fprintf(os.Stderr, "\nOnly one -pattern command is permitted\n")
		os.Exit(1)
	}

	// look for improperly formed commands
	ValidateCommands(head)

	if noElement {
		fmt.Fprintf(os.Stderr, "\nNo -element statement in argument list\n")
		os.Exit(1)
	}

	// convert command strings to array of operations for faster processing
	ResolveOperations(head)

	if len(head.Subtasks) != 1 {
		return nil
	}

	// skip past empty placeholder
	return head.Subtasks[0]
}

// COLLECT AND FORMAT REQUESTED XML VALUES

// parse attributes only if attribute value is requested in element statement
func ParseAttributes(attrib string) []string {

	if attrib == "" {
		return nil
	}

	// count equal signs
	num := 0
	for i := 0; i < len(attrib); i++ {
		if attrib[i] == '=' {
			num += 2
		}
	}
	if num < 1 {
		return nil
	}

	// allocate array of proper size
	arry := make([]string, num)
	if arry == nil {
		return nil
	}

	start := 0
	idx := 0
	itm := 0

	// place tag and value in successive array slots
	for idx < len(attrib) && itm < num {
		ch := attrib[idx]
		if ch == '=' {
			// =
			arry[itm] = attrib[start:idx]
			itm++
			// skip past equal sign and leading double quote
			idx += 2
			start = idx
		} else if ch == '"' {
			// "
			arry[itm] = attrib[start:idx]
			itm++
			// skip past trailing double quote and (possible) space
			idx += 2
			start = idx
		} else {
			idx++
		}
	}

	return arry
}

// return matched element values to callback
func ExploreElements(curr *Node, prnt, match, attrib string, level int, proc func(string, int)) {

	if curr == nil || proc == nil {
		return
	}

	if (curr.Name == match || (match == "" && attrib != "")) && (prnt == "" || curr.Parent == prnt) {

		if attrib != "" {
			if curr.Attributes != "" && curr.Attribs == nil {
				// parse attributes on-the-fly if queried
				curr.Attribs = ParseAttributes(curr.Attributes)
			}
			for i := 0; i < len(curr.Attribs)-1; i += 2 {
				// attributes now parsed into array as [ tag, value, tag, value, tag, value, ... ]
				if curr.Attribs[i] == attrib {
					proc(curr.Attribs[i+1], level)
					return
				}
			}

		} else if curr.Contents != "" {

			proc(curr.Contents, level)
			return

		} else if curr.Children != nil {

			// for XML container object, send empty string to callback to increment count
			proc("", level)
			// and continue exploring
		}
	}

	for chld := curr.Children; chld != nil; chld = chld.Next {
		ExploreElements(chld, prnt, match, attrib, level+1, proc)
	}
}

const (
	_ = iota
	COMPACT
	FLUSH
	INDENT
)

// -element "*" prints XML subtree
func PrintSubtree(node *Node, style int, printAttrs bool, proc func(string)) {

	if node == nil || proc == nil {
		return
	}

	// array to speed up indentation
	indentSpaces := []string{
		"",
		"  ",
		"    ",
		"      ",
		"        ",
		"          ",
		"            ",
		"              ",
		"                ",
		"                  ",
	}

	var doSubtree func(*Node, int)

	// doSubtree must be defined in var section in order to call itself recursively
	doSubtree = func(curr *Node, depth int) {

		// suppress if it would be an empty self-closing tag
		if !IsNotJustWhitespace(curr.Attributes) && curr.Contents == "" && curr.Children == nil {
			return
		}

		if style == INDENT {
			i := depth
			for i > 9 {
				proc("                    ")
				i -= 10
			}
			proc(indentSpaces[i])
		}

		proc("<")
		proc(curr.Name)

		if printAttrs {
			attr := strings.TrimSpace(curr.Attributes)
			attr = CompressRunsOfSpaces(attr)
			if attr != "" {
				proc(" ")
				proc(attr)
			}
		}

		// see if suitable for for self-closing tag
		if curr.Contents == "" && curr.Children == nil {
			proc("/>")
			if style != COMPACT {
				proc("\n")
			}
			return
		}

		proc(">")

		if curr.Contents != "" {

			str := html.EscapeString(curr.Contents)
			proc(str)

		} else {

			if style != COMPACT {
				proc("\n")
			}

			for chld := curr.Children; chld != nil; chld = chld.Next {
				doSubtree(chld, depth+1)
			}

			if style == INDENT {
				i := depth
				for i > 9 {
					proc("                    ")
					i -= 10
				}
				proc(indentSpaces[i])
			}
		}

		proc("<")
		proc("/")
		proc(curr.Name)
		proc(">")

		if style != COMPACT {
			proc("\n")
		}
	}

	doSubtree(node, 1)
}

func ProcessElement(curr *Node, item string, status, index, level int, variables map[string]string, acc func(string)) {

	if curr == nil || acc == nil {
		return
	}

	style := COMPACT
	printAttrs := true

	// check for special character at beginning of name
	if len(item) > 1 {
		switch item[0] {
		case '&':
			if IsAllCapsOrDigits(item[1:]) {
				status = VARIABLE
				item = item[1:]
			} else {
				fmt.Fprintf(os.Stderr, "\nUnrecognized variable '%s'\n", item)
				os.Exit(1)
			}
		case '#':
			status = COUNT
			item = item[1:]
		case '%':
			status = LENGTH
			item = item[1:]
		case '^':
			status = DEPTH
			item = item[1:]
		case '*':
			status = STAR
			item = item[1:]
			for _, ch := range item {
				if ch == '*' {
					style++
				} else if ch == '@' {
					printAttrs = false
				}
			}
			if style > INDENT {
				style = INDENT
			}
		default:
		}
	} else if item == "*" {
		status = STAR
	} else if item == "+" {
		// -element "+" prints index of current XML object
		val := strconv.Itoa(index)
		acc(val)
		return
	}

	match, val := SplitInTwoAt(item, ":", LEFT)
	prnt, match := SplitInTwoAt(match, "/", RIGHT)
	match, attrib := SplitInTwoAt(match, "@", LEFT)

	// element selection should not have colon followed by value
	if val != "" {
		fmt.Fprintf(os.Stderr, "\nelement-colon-value construct '%s' is inappropriate for item selection\n", item)
		os.Exit(1)
	}

	switch status {
	case ELEMENT, ENCODE:
		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				if str != "" {
					acc(str)
				}
			})
	case VARIABLE:
		// use value of stored variable
		val, ok := variables[match]
		if ok {
			acc(val)
		}
	case COUNT:
		count := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				count++
			})

		// number of element objects
		val := strconv.Itoa(count)
		acc(val)
	case LENGTH:
		length := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				length += len(str)
			})

		// length of element strings
		val := strconv.Itoa(length)
		acc(val)
	case DEPTH:
		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				// depth of each element in scope
				val := strconv.Itoa(lvl)
				acc(val)
			})
	case FIRST:
		single := ""

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				if single == "" {
					single = str
				}
			})

		if single != "" {
			acc(single)
		}
	case LAST:
		single := ""

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				single = str
			})

		if single != "" {
			acc(single)
		}
	case SUM:
		sum := 0
		ok := false

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				value, err := strconv.Atoi(str)
				if err == nil {
					sum += value
					ok = true
				}
			})

		if ok {
			// sum of element values
			val := strconv.Itoa(sum)
			acc(val)
		}
	case MIN:
		min := 0
		ok := false

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				value, err := strconv.Atoi(str)
				if err == nil {
					if !ok || value < min {
						min = value
					}
					ok = true
				}
			})

		if ok {
			// minimum of element values
			val := strconv.Itoa(min)
			acc(val)
		}
	case MAX:
		max := 0
		ok := false

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				value, err := strconv.Atoi(str)
				if err == nil {
					if !ok || value > max {
						max = value
					}
					ok = true
				}
			})

		if ok {
			// maximum of element values
			val := strconv.Itoa(max)
			acc(val)
		}
	case STAR:
		// -element "*" prints current XML subtree on a single line
		var buffer bytes.Buffer

		PrintSubtree(curr, style, printAttrs,
			func(str string) {
				if str != "" {
					buffer.WriteString(str)
				}
			})

		txt := buffer.String()
		if txt != "" {
			acc(txt)
		}
	default:
	}
}

func ProcessClause(curr *Node, str, prev, pfx, sfx, sep string, status, index, level int, variables map[string]string) (string, bool) {

	ok := false
	num := 0

	// format results in buffer
	var buffer bytes.Buffer

	buffer.WriteString(prev)
	buffer.WriteString(pfx)
	between := ""

	// element names combined with commas are treated as a prefix-separator-suffix group
	comma := strings.Split(str, ",")
	for _, item := range comma {

		ProcessElement(curr, item, status, index, level, variables,
			func(str string) {
				if str != "" {
					switch status {
					case ELEMENT, FIRST, LAST:
						buffer.WriteString(between)
						buffer.WriteString(str)
						between = sep
					case ENCODE:
						buffer.WriteString(between)
						str = html.EscapeString(str)
						buffer.WriteString(str)
						between = sep
					case SUM, MIN, MAX:
						value, err := strconv.Atoi(str)
						if err != nil {
							break
						}
						switch status {
						case SUM:
							num += value
						case MIN:
							if !ok || value < num {
								num = value
							}
						case MAX:
							if !ok || value > num {
								num = value
							}
						default:
						}
					default:
					}

					ok = true
				}
			})
	}

	if status == SUM || status == MIN || status == MAX {
		if ok {
			val := strconv.Itoa(num)
			buffer.WriteString(val)
		}
	}

	buffer.WriteString(sfx)

	if !ok {
		return "", false
	}

	txt := buffer.String()

	return txt, true
}

func ProcessInstructions(commands []*Operation, curr *Node, tab, ret string, index, level int, variables map[string]string, accum func(string)) (string, string) {

	if accum == nil {
		return tab, ret
	}

	sep := "\t"
	pfx := ""
	sfx := ""

	col := "\t"
	lin := "\n"

	varname := ""

	// process commands
	for _, op := range commands {
		str := op.Value
		switch op.Type {
		case ELEMENT, FIRST, LAST, ENCODE, SUM, MIN, MAX:
			txt, ok := ProcessClause(curr, str, tab, pfx, sfx, sep, op.Type, index, level, variables)
			if ok {
				tab = col
				ret = lin
				accum(txt)
			}
		case TAB:
			col = str
		case RET:
			lin = str
		case PFX:
			pfx = str
		case SFX:
			sfx = str
		case SEP:
			sep = str
		case LBL:
			lbl := str
			accum(tab)
			accum(lbl)
			tab = col
			ret = lin
		case PFC:
			// preface clears previous tab and sets prefix in one command
			pfx = str
			fallthrough
		case CLR:
			// clear previous tab after the fact
			tab = ""
		case RST:
			pfx = ""
			sfx = ""
			sep = "\t"
		case VARIABLE:
			varname = str
		case VALUE:
			len := len(str)
			if len > 1 && str[0] == '(' && str[len-1] == ')' {
				// set variable from literal text inside parentheses, e.g., -COM "(, )"
				variables[varname] = str[1 : len-1]
				// -match "&VARIABLE" will succeed if set to blank with empty parentheses "()"
			} else if str == "" {
				// -match "&VARIABLE" will fail if initialized with empty string ""
				delete(variables, varname)
			} else {
				txt, ok := ProcessClause(curr, str, "", pfx, sfx, sep, ELEMENT, index, level, variables)
				if ok {
					variables[varname] = txt
				}
			}
			varname = ""
		default:
		}
	}

	return tab, ret
}

// -MATCH AND -AVOID CONDITIONAL TESTS

func HasElement(curr *Node, str string, level int, variables map[string]string) bool {

	match, val := SplitInTwoAt(str, ":", LEFT)
	prnt, match := SplitInTwoAt(match, "/", RIGHT)
	match, attrib := SplitInTwoAt(match, "@", LEFT)

	found := false

	isVariable := false

	// skip past pound, percent, or caret character at beginning of name
	if len(match) > 1 {
		switch match[0] {
		case '&':
			if IsAllCapsOrDigits(match[1:]) {
				isVariable = true
				match = match[1:]
			} else {
				fmt.Fprintf(os.Stderr, "\nUnrecognized variable '%s'\n", match)
				os.Exit(1)
			}
		case '#', '%', '^':
			match = match[1:]
		default:
		}
	} else if match == "+" {
		return true
	}

	if isVariable {

		// use value of stored variable
		str, ok := variables[match]
		if ok {
			if val == "" || val == str {
				found = true
			}
		}

	} else {

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				if val == "" || val == str {
					found = true
				}
			})
	}

	return found
}

func GetValue(curr *Node, str string, index, level int, variables map[string]string) (int, bool) {

	// first character may be backslash protecting minus sign (undocumented)
	if len(str) > 1 && str[0] == '\\' {
		str = str[1:]
	}

	if str != "+" {
		// check for numeric argument
		number, err := strconv.Atoi(str)
		if err == nil {
			return number, true
		}
	}

	result := 0
	found := false

	status := ELEMENT

	// check for pound, percent, or caret character at beginning of name
	if len(str) > 1 {
		switch str[0] {
		case '&':
			if IsAllCapsOrDigits(str[1:]) {
				status = VARIABLE
				str = str[1:]
			} else {
				fmt.Fprintf(os.Stderr, "\nUnrecognized variable '%s'\n", str)
				os.Exit(1)
			}
		case '#':
			status = COUNT
			str = str[1:]
		case '%':
			status = LENGTH
			str = str[1:]
		case '^':
			status = DEPTH
			str = str[1:]
		default:
		}
	} else if str == "+" {
		// index of explored parent object
		return index, true
	}

	match, val := SplitInTwoAt(str, ":", LEFT)
	prnt, match := SplitInTwoAt(match, "/", RIGHT)
	match, attrib := SplitInTwoAt(match, "@", LEFT)

	// -match argument for numeric comparison should not have colon followed by value
	if val != "" {
		fmt.Fprintf(os.Stderr, "\nelement-colon-value construct '%s' is inappropriate for numeric comparison\n", str)
		os.Exit(1)
	}

	switch status {
	case ELEMENT:
		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				value, err := strconv.Atoi(str)
				if err == nil {
					result = value
					found = true
				}
			})
	case VARIABLE:
		// use value of stored variable
		str, ok := variables[match]
		if ok {
			value, err := strconv.Atoi(str)
			if err == nil {
				result = value
				found = true
			}
		}
	case COUNT:
		count := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				count++
				found = true
			})

		// number of element objects
		result = count
	case LENGTH:
		length := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				length += len(str)
				found = true
			})

		// length of element strings
		result = length
	case DEPTH:
		depth := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				depth = lvl
				found = true
			})

		// depth of last element in scope
		result = depth
	default:
	}

	return result, found
}

func InRange(curr *Node, str, value string, status, index, level int, variables map[string]string) bool {

	// first argument must be an element/count/length/variable (numeric value will fail original -match)
	x, okx := GetValue(curr, str, index, level, variables)

	// second argument (after -lt, -ge, etc.) may be numeric, but can also be an element expression
	y, oky := GetValue(curr, value, index, level, variables)

	// both arguments must resolve to integers
	if !okx || !oky {
		return false
	}

	switch status {
	case GT:
		if x > y {
			return true
		}
	case GE:
		if x >= y {
			return true
		}
	case LT:
		if x < y {
			return true
		}
	case LE:
		if x <= y {
			return true
		}
	case EQ:
		if x == y {
			return true
		}
	case NE:
		if x != y {
			return true
		}
	default:
	}

	return false
}

func GetString(curr *Node, str string, index, level int, variables map[string]string) (string, bool) {

	result := ""
	found := false

	status := ELEMENT

	// check for pound, percent, or caret character at beginning of name
	if len(str) > 1 {
		switch str[0] {
		case '&':
			if IsAllCapsOrDigits(str[1:]) {
				status = VARIABLE
				str = str[1:]
			} else {
				fmt.Fprintf(os.Stderr, "\nUnrecognized variable '%s'\n", str)
				os.Exit(1)
			}
		case '#':
			status = COUNT
			str = str[1:]
		case '%':
			status = LENGTH
			str = str[1:]
		case '^':
			status = DEPTH
			str = str[1:]
		default:
		}
	} else if str == "+" {
		// index of explored parent object
		return strconv.Itoa(index), true
	}

	match, val := SplitInTwoAt(str, ":", LEFT)
	prnt, match := SplitInTwoAt(match, "/", RIGHT)
	match, attrib := SplitInTwoAt(match, "@", LEFT)

	// -match argument for substring comparison should not have colon followed by value
	if val != "" {
		fmt.Fprintf(os.Stderr, "\nelement-colon-value construct '%s' is inappropriate for substring comparison\n", str)
		os.Exit(1)
	}

	switch status {
	case ELEMENT:
		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				if str != "" {
					result = str
					found = true
				}
			})
	case VARIABLE:
		// use value of stored variable
		str, ok := variables[match]
		if ok {
			result = str
			found = true
		}
	case COUNT:
		count := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				count++
				found = true
			})

		// number of element objects
		result = strconv.Itoa(count)
	case LENGTH:
		length := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				length += len(str)
				found = true
			})

		// length of element strings
		result = strconv.Itoa(length)
	case DEPTH:
		depth := 0

		ExploreElements(curr, prnt, match, attrib, level,
			func(str string, lvl int) {
				depth = lvl
				found = true
			})

		// depth of last element in scope
		result = strconv.Itoa(depth)
	default:
	}

	return result, found
}

func InString(curr *Node, str, substr string, status, index, level int, variables map[string]string) bool {

	if substr == "" {
		return false
	}

	// first character may be backslash protecting minus sign (undocumented)
	if len(substr) > 1 && substr[0] == '\\' {
		substr = substr[1:]
	}

	// first argument must be an element/count/length/variable)
	val, ok := GetString(curr, str, index, level, variables)
	if !ok {
		return false
	}

	switch status {
	case EQUALS:
		if strings.ToUpper(val) == strings.ToUpper(substr) {
			return true
		}
	case CONTAINS:
		if strings.Contains(strings.ToUpper(val), strings.ToUpper(substr)) {
			return true
		}
	case STARTSWITH:
		if strings.HasPrefix(strings.ToUpper(val), strings.ToUpper(substr)) {
			return true
		}
	case ENDSWITH:
		if strings.HasSuffix(strings.ToUpper(val), strings.ToUpper(substr)) {
			return true
		}
	default:
	}

	return false
}

func ConditionsAreSatisfied(conditions []*Operation, curr *Node, index, level int, variables map[string]string) bool {

	required := 0
	observed := 0
	forbidden := 0
	isMatch := false
	isAvoid := false

	// previous element name needs to be remembered for subsequent numeric or substring tests
	prev := ""

	// test conditional arguments
	for _, op := range conditions {
		str := op.Value
		switch op.Type {
		// -match tests for presence of element (or element with specific value)
		case MATCH:
			// checking for failure here allows for multiple -match [ -and / -or ] clauses
			if isMatch && observed < required {
				return false
			}
			if isAvoid && forbidden > 0 {
				return false
			}
			required = 0
			observed = 0
			forbidden = 0
			isMatch = true
			isAvoid = false
			// continue on to next two cases
			fallthrough
		case AND:
			required++
			// continue on to next case
			fallthrough
		case OR:
			prev = str
			if HasElement(curr, str, level, variables) {
				observed++
				// record presence of forbidden element if in -avoid clause
				forbidden++
			}
		// -avoid tests for absence of element (or element with specific value)
		case AVOID:
			if isMatch && observed < required {
				return false
			}
			if isAvoid && forbidden > 0 {
				return false
			}
			required = 0
			observed = 0
			forbidden = 0
			isMatch = false
			isAvoid = true
			prev = str
			if HasElement(curr, str, level, variables) {
				// record presence of forbidden element, can be negated by failure of subsequent substring or numeric constraint command
				forbidden++
			}
		case EQUALS, CONTAINS, STARTSWITH, ENDSWITH:
			// substring test on element values for -match
			required++
			if InString(curr, prev, str, op.Type, index, level, variables) {
				observed++
				forbidden++
			} else {
				// neutralize presence of element in -avoid statement
				forbidden--
			}
		case GT, GE, LT, LE, EQ, NE:
			// numeric tests on element values for -match
			required++
			if InRange(curr, prev, str, op.Type, index, level, variables) {
				observed++
				forbidden++
			} else {
				// failure of conditional suppresses -avoid element match
				forbidden--
			}
		default:
		}
	}

	if isMatch && observed < required {
		return false
	}
	if isAvoid && forbidden > 0 {
		return false
	}

	return true
}

// RECURSIVELY PROCESS EXPLORATION COMMANDS AND XML DATA STRUCTURE

func ExploreNodes(curr *Node, prnt, match string, index, level int, proc func(*Node, int, int)) int {

	if curr == nil || proc == nil {
		return index
	}

	// match is "*" for heterogeneous data constructs, e.g., -group PubmedArticleSet/*
	if (curr.Name == match || match == "*") && (prnt == "" || curr.Parent == prnt) {
		proc(curr, index, level)
		return index + 1
	}

	// clearing prnt "*" now allows exploration within recursive data, e.g., -pattern Taxon -block */Taxon
	if prnt == "*" {
		prnt = ""
	}

	// explore child nodes
	for chld := curr.Children; chld != nil; chld = chld.Next {
		index = ExploreNodes(chld, prnt, match, index, level+1, proc)
	}

	return index
}

func ProcessCommands(cmds *Block, curr *Node, tab, ret string, index, level int, variables map[string]string, accum func(string)) (string, string) {

	if accum == nil {
		return tab, ret
	}

	if len(cmds.Commands) > 0 {
		return ProcessInstructions(cmds.Commands, curr, tab, ret, index, level, variables, accum)
	}

	if cmds.Visit == "" {
		return tab, ret
	}

	// closure passes local variables to callback
	processNode := func(node *Node, idx, lvl int) {

		// apply -match or -avoid tests
		if ConditionsAreSatisfied(cmds.Conditions, node, idx, lvl, variables) {

			// process sub commands on child node
			for _, sub := range cmds.Subtasks {
				tab, ret = ProcessCommands(sub, node, tab, ret, idx, lvl, variables, accum)
			}
		}
	}

	prnt, match := SplitInTwoAt(cmds.Visit, "/", RIGHT)

	// apply -position test
	if cmds.Position == "" {

		ExploreNodes(curr, prnt, match, 1, level, processNode)

	} else if cmds.Position == "first" {

		var single *Node
		lev := 0
		ind := 0

		ExploreNodes(curr, prnt, match, 1, level,
			func(node *Node, idx, lvl int) {
				if single == nil {
					single = node
					ind = idx
					lev = lvl
				}
			})

		processNode(single, ind, lev)

	} else if cmds.Position == "last" {

		var single *Node
		lev := 0
		ind := 0

		ExploreNodes(curr, prnt, match, 1, level,
			func(node *Node, idx, lvl int) {
				single = node
				ind = idx
				lev = lvl
			})

		processNode(single, ind, lev)

	} else {

		// use numeric position
		number, err := strconv.Atoi(cmds.Position)
		if err == nil {

			pos := 0

			ExploreNodes(curr, prnt, match, 1, level,
				func(node *Node, idx, lvl int) {
					pos++
					if pos == number {
						processNode(node, pos, lvl)
					}
				})

		} else {

			fmt.Fprintf(os.Stderr, "\nUnrecognized position '%s'\n", cmds.Position)
			os.Exit(1)
		}
	}

	return tab, ret
}

// PROCESS ONE XML COMPONENT RECORD

func ParseXml(start, attrib, parent string, rdr *Tokenizer) (*Node, bool) {

	if rdr == nil {
		return nil, false
	}

	ok := true
	node := &Node{Name: start, Parent: parent, Attributes: attrib}

	var lastNode *Node

	for {
		// read next token
		cat, name, attributes, contents := NextToken(rdr)

		switch cat {
		case IS_START_TAG:
			// read sub tree
			obj, ok := ParseXml(name, attributes, node.Name, rdr)
			if !ok {
				break
			}
			// adding next child to end of linked list gives better performance than appending to slice of nodes
			if node.Children == nil {
				node.Children = obj
			}
			if lastNode != nil {
				lastNode.Next = obj
			}
			lastNode = obj
		case IS_END_TAG:
			// pop out of recursive call
			return node, ok
		case IS_CONTENT:
			node.Contents = contents
		case IS_CLOSED:
			return nil, false
		default:
		}
	}
}

func ProcessQuery(start, attrib, parent string, cmds *Block, rdr *Tokenizer) (string, bool) {

	if cmds == nil || rdr == nil {
		return "", false
	}

	// exit from function will collect garbage of node structure for current XML object
	pat, ok := ParseXml(start, attrib, parent, rdr)

	if !ok {
		return "", false
	}

	// exit from function will also free map of recorded variables for current -pattern
	variables := make(map[string]string)

	var buffer bytes.Buffer

	// start processing at top of command tree and top of XML subregion selected by -pattern
	_, ret := ProcessCommands(cmds, pat, "", "", 1, 1, variables,
		func(str string) {
			if str != "" {
				buffer.WriteString(str)
			}
		})

	if ret != "" {
		buffer.WriteString(ret)
	}

	txt := buffer.String()

	// return consolidated result string
	return txt, true
}

// CONCURRENT PRODUCER AND CONSUMER GOROUTINES PARSE AND PROCESS PARTITIONED XML OBJECTS

type Extract struct {
	Index int
	Text  string
}

// producer sends partitioned XML strings through channel
func XmlProducer(pat, star string, rdr *XmlReader, out chan Extract) {

	// close channel when all records have been processed, so consumers can range over channel
	defer close(out)

	// partition all input by pattern and send XML substring to available consumer through channel
	PartitionPattern(pat, star, rdr,
		func(rec int, str string) {
			out <- Extract{rec, str}
		})
}

func CreateProducer(pat, star string, chan_depth int, rdr *XmlReader) chan Extract {

	out := make(chan Extract, chan_depth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create producer channel\n")
		os.Exit(1)
	}

	// launch single producer goroutine
	go XmlProducer(pat, star, rdr, out)

	return out
}

// consumer reads partitioned XML from channel and calls parser for processing
func ScanAndProcess(text, parent string, cmds *Block) string {

	// create XML parser
	rdr := NewTokenizer(nil, text)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create XML Tokenizer\n")
		os.Exit(1)
	}

	str := ""

	// scan xml tags, process when first start tag is found
	for {
		cat, name, attributes, _ := NextToken(rdr)
		if cat == IS_CLOSED {
			return ""
		}
		if cat != IS_START_TAG {
			continue
		}

		// read and process one -pattern object at a time
		str, _ = ProcessQuery(name, attributes, parent, cmds, rdr)

		break
	}

	return str
}

func XmlConsumer(cmds *Block, parent string, wg *sync.WaitGroup, inp chan Extract, out chan Extract) {

	// report when this consumer has no more records to process
	defer wg.Done()

	// read partitioned XML from producer channel
	for ext := range inp {

		idx := ext.Index
		text := ext.Text

		if text == "" {
			// should never see empty input data
			out <- Extract{idx, text}
			continue
		}

		str := ScanAndProcess(text, parent, cmds)

		// send even if empty to get all record counts for reordering
		out <- Extract{idx, str}
	}
}

func CreateConsumers(cmds *Block, parent string, num_servers, chan_depth int, inp chan Extract) chan Extract {

	out := make(chan Extract, chan_depth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create consumer channel\n")
		os.Exit(1)
	}

	var wg sync.WaitGroup

	// launch multiple consumer goroutines
	for i := 0; i < num_servers; i++ {
		wg.Add(1)
		go XmlConsumer(cmds, parent, &wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all consumers are done, then close single output channel, so unshuffler can range over channel
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// unshuffler restores output of multiple consumers to original record order
type ExtractHeap []Extract

// methods that satisfy heap.Interface
func (h ExtractHeap) Len() int {
	return len(h)
}
func (h ExtractHeap) Less(i, j int) bool {
	return h[i].Index < h[j].Index
}
func (h ExtractHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}
func (h *ExtractHeap) Push(x interface{}) {
	*h = append(*h, x.(Extract))
}
func (h *ExtractHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func Unshuffler(heap_size int, inp chan Extract, out chan Extract) {

	// close channel when all records have been processed
	defer close(out)

	// initialize empty heap
	hp := &ExtractHeap{}
	heap.Init(hp)

	// index of next desired result
	next := 1

	delay := 0

	for ext := range inp {

		// push result onto heap
		heap.Push(hp, ext)

		// read several values before checking to see if next record to print has been processed
		if delay < heap_size {
			delay++
			continue
		}

		delay = 0

		for hp.Len() > 0 {
			// remove lowest item from heap, use interface type assertion
			curr := heap.Pop(hp).(Extract)

			if curr.Index == next {
				// if this is the desired item, send to output channel
				out <- curr
				// increment index
				next++
				// and keep checking heap to see if next result is already available
			} else {
				// otherwise push back onto heap
				heap.Push(hp, curr)
				// and go back to waiting on input channel
				break
			}
		}
	}

	// send remainder of heap to output channel
	for hp.Len() > 0 {
		out <- heap.Pop(hp).(Extract)
	}
}

func CreateUnshuffler(chan_depth, heap_size int, inp chan Extract) chan Extract {

	out := make(chan Extract, chan_depth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create unshuffler channel\n")
		os.Exit(1)
	}

	// launch single unshuffler goroutine
	go Unshuffler(heap_size, inp, out)

	return out
}

// ORIGINAL STREAMER/PARSER CODE

// function ProcessXml is obsolete, but is retained in case some unusual situation is able to cause a goroutine deadlock
func ProcessXml(cmds *Block, topPat, star string, blk *XmlReader) int {

	if cmds == nil || topPat == "" || blk == nil {
		return 0
	}

	count := 0

	// create XML parser
	rdr := NewTokenizer(blk, "")
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create XML Tokenizer\n")
		os.Exit(1)
	}

	// scan xml tags, process as soon as each top-level pattern is read
	for {
		cat, name, attributes, _ := NextToken(rdr)
		if cat == IS_CLOSED {
			return count
		}
		if cat != IS_START_TAG {
			continue
		}
		if name != topPat {
			continue
		}

		if star == "" {
			// read and process one -pattern object at a time
			str, ok := ProcessQuery(name, attributes, "", cmds, rdr)
			if !ok {
				return count
			}
			count++
			if str != "" {
				fmt.Printf("%s", str)
			}
		} else {
			// read and process heterogeneous objects immediately below -pattern parent
			for {
				cat, name, attributes, _ := NextToken(rdr)
				if cat == IS_CLOSED {
					return count
				}
				if cat != IS_START_TAG {
					// break to scan for next -pattern parent if XML downloaded in chunks
					break
				}
				str, ok := ProcessQuery(name, attributes, topPat, cmds, rdr)
				if !ok {
					return count
				}
				count++
				if str != "" {
					fmt.Printf("%s", str)
				}
			}
		}
	}
}

// MAIN FUNCTION

// e.g., xtract -pattern PubmedArticle -element MedlineCitation/PMID -block Author -element Initials,LastName

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nNo command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// DOCUMENTATION COMMANDS

	inSwitch := true

	switch args[0] {
	case "-version":
		fmt.Printf("%s\n", xtract_version)
	case "-help":
		fmt.Printf("xtract %s\n%s\n", xtract_version, xtract_help)
	case "-extras", "-extra":
		fmt.Printf("xtract %s\n%s\n", xtract_version, xtract_extras)
	case "-examples", "-example":
		fmt.Printf("xtract %s\n%s\n", xtract_version, xtract_examples)
	case "-scripts", "-script":
		fmt.Printf("xtract %s\n%s\n", xtract_version, xtract_scripts)
	case "-internal", "-internals":
		fmt.Printf("xtract %s\n%s\n", xtract_version, xtract_internal)
	case "-keys":
		fmt.Printf("%s\n", keyboard_shortcuts)
	case "-unix":
		fmt.Printf("%s\n", unix_commands)
	default:
		// if not any of the documentation commands, keep going
		inSwitch = false
	}

	if inSwitch {
		return
	}

	// CONCURRENCY, CLEANUP, AND DEBUGGING FLAGS

	ncpu := runtime.NumCPU()

	// access to earlier parsing mechanism just in case
	legacy := false

	// concurrent performance tuning defaults
	num_procs := 8
	num_servers := 32
	chan_depth := 8
	heap_size := 16

	// adjust defaults for large multi-processor computer
	if ncpu > 16 {
		num_procs = 16
		num_servers = 64
	}

	// limit processors to actual hardware
	if num_procs > ncpu {
		num_procs = ncpu
	}

	// XML data cleanup
	doCompress := false
	doCleanup := false

	// debugging
	timer := false
	debug := false
	empty := false
	index := false

	// function to get numeric value
	getNumericArg := func(name string, zer, max int) int {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\n%s is missing\n", name)
			os.Exit(1)
		}
		value, err := strconv.Atoi(args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "\n%s (%s) is not an integer\n", name, args[1])
			os.Exit(1)
		}
		// skip past first of two arguments
		args = args[1:]

		// special case for argument value of 0
		if value < 1 {
			return zer
		}
		// limit value to specified maximum
		if value > max {
			return max
		}
		return value
	}

	// get concurrency, cleanup, and debugging flags in any order
	for {

		inSwitch = true

		switch args[0] {
		case "-legacy":
			legacy = true
		case "-proc":
			num_procs = getNumericArg("Number of processors", ncpu, ncpu)
		case "-serv":
			num_servers = getNumericArg("Number of servers", ncpu, 128)
		case "-chan":
			chan_depth = getNumericArg("Depth of channels", ncpu, 32)
		case "-heap":
			heap_size = getNumericArg("Size of heap", 16, 32)
		case "-compress":
			doCompress = true
		case "-cleanup":
			doCleanup = true
		case "-timer":
			timer = true
		case "-debug":
			debug = true
			timer = true
		case "-empty":
			empty = true
		case "-index":
			index = true
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past argument
		args = args[1:]

		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "\nInsufficient command-line arguments supplied to xtract\n")
			os.Exit(1)
		}
	}

	// allow simultaneous go routine threads
	runtime.GOMAXPROCS(num_procs)

	if debug {
		fmt.Fprintf(os.Stderr, "\nXtract\n\n")
		fmt.Fprintf(os.Stderr, "  CPU:  %d\n  Proc: %d\n  Serv: %d\n  Chan: %d\n  Heap: %d\n\n", ncpu, num_procs, num_servers, chan_depth, heap_size)
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nInsufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// SPECIAL FORMATTING COMMANDS

	inSwitch = true

	switch args[0] {
	case "-outline":
		ProcessOutline(doCompress, doCleanup)
	case "-synopsis":
		ProcessSynopsis(doCompress, doCleanup)
	case "-verify", "-validate":
		ProcessVerify(doCompress, doCleanup)
	default:
		// if not any of the formatting commands, keep going
		inSwitch = false
	}

	if inSwitch {
		return
	}

	// INITIALIZE PROCESS TIMER AND RECORD COUNT

	startTime := time.Now()
	recordCount := 0
	byteCount := 0

	// function to print processing rate and program duration
	printDuration := func(name string) {
		endTime := time.Now()
		duration := endTime.Sub(startTime)
		seconds := float64(duration.Nanoseconds()) / 1e9

		fmt.Fprintf(os.Stderr, "\nXtract processed %d %s in %.3f seconds", recordCount, name, seconds)
		if seconds >= 0.001 && recordCount > 0 {
			rate := int(0.5 + float64(recordCount)/seconds)
			fmt.Fprintf(os.Stderr, " (%d %s/second", rate, name)
			if byteCount > 0 {
				rate := int(0.5 + float64(byteCount)/seconds)
				fmt.Fprintf(os.Stderr, ", %d bytes/second", rate)
			}
			fmt.Fprintf(os.Stderr, ")")
		}
		fmt.Fprintf(os.Stderr, "\n\n")
	}

	// PERFORMANCE TIMING COMMANDS

	// -chunk tests speed of NextBlock streaming
	if args[0] == "-chunk" {
		in := bufio.NewReader(os.Stdin)
		blk := NewXmlReader(in, doCompress, doCleanup)
		if blk == nil {
			return
		}

		for {
			str := NextBlock(blk)
			if str == "" {
				break
			}
			recordCount++
			byteCount += len(str)
		}

		printDuration("blocks")

		return
	}

	// -split tests speed of PartitionPattern streaming, requires pattern
	if args[0] == "-split" {
		if len(args) > 1 {
			if args[1] == "-pattern" || args[1] == "-Pattern" {
				// skip past -split if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nPattern missing after -split command\n")
			os.Exit(1)
		}
		pat := args[1]

		in := bufio.NewReader(os.Stdin)
		blk := NewXmlReader(in, doCompress, doCleanup)
		if blk == nil {
			return
		}

		PartitionPattern(pat, "", blk,
			func(rec int, str string) {
				recordCount++
				byteCount += len(str)
			})

		printDuration("patterns")

		return
	}

	// -token tests speed of NextToken streaming, requires pattern
	if args[0] == "-token" {
		if len(args) > 1 {
			if args[1] == "-pattern" || args[1] == "-Pattern" {
				// skip past -token if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nPattern missing after -token command\n")
			os.Exit(1)
		}
		pat := args[1]

		in := bufio.NewReader(os.Stdin)
		blk := NewXmlReader(in, doCompress, doCleanup)
		if blk == nil {
			return
		}

		PartitionPattern(pat, "", blk,
			func(rec int, str string) {
				rdr := NewTokenizer(nil, str)
				if rdr == nil {
					return
				}
				for {
					cat, name, attributes, contents := NextToken(rdr)
					if cat == IS_CLOSED {
						break
					}
					recordCount++
					byteCount += len(name) + len(attributes) + len(contents)
				}
			})

		printDuration("tokens")

		return
	}

	// -parse tests speed of ProcessXML streaming, requires pattern
	if args[0] == "-parse" {
		if len(args) > 1 {
			if args[1] == "-pattern" || args[1] == "-Pattern" {
				// skip past -parse if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nPattern missing after -parse command\n")
			os.Exit(1)
		}
		pat := args[1]

		in := bufio.NewReader(os.Stdin)
		blk := NewXmlReader(in, doCompress, doCleanup)
		if blk == nil {
			return
		}

		PartitionPattern(pat, "", blk,
			func(rec int, str string) {
				rdr := NewTokenizer(nil, str)
				if rdr == nil {
					return
				}
				for {
					cat, name, attributes, _ := NextToken(rdr)
					if cat == IS_CLOSED {
						break
					}
					if cat != IS_START_TAG {
						continue
					}
					ParseXml(name, attributes, "", rdr)
					recordCount++
					byteCount += len(str)
				}
			})

		printDuration("objects")

		return
	}

	// SEQUENCE RECORD EXTRACTION COMMAND GENERATOR

	// -insd simplifies extraction of INSDSeq qualifiers
	if args[0] == "-insd" {

		args = args[1:]

		fi, _ := os.Stdin.Stat()
		isPipe := bool((fi.Mode() & os.ModeCharDevice) == 0)

		insd := ProcessINSD(args, isPipe)

		if !isPipe {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range insd {
				fmt.Printf(" %s", str)
			}
			// fmt.Printf("| \\\n")
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = insd
	}

	// CITATION MATCHER EXTRACTION COMMAND GENERATOR

	// -hydra filters HydraResponse output by relevance score (undocumented)
	if args[0] == "-hydra" {

		fi, _ := os.Stdin.Stat()
		isPipe := bool((fi.Mode() & os.ModeCharDevice) == 0)

		hydra := ProcessHydra(isPipe)

		if !isPipe {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range hydra {
				fmt.Printf(" %s", str)
			}
			// fmt.Printf("| \\\n")
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = hydra
	}

	// XML DATA FORMATTING/COMPRESSION

	// -format can take an optional parent pattern and put each XML sub-object on its own line, for fastest subsequent processing (undocumented)
	if args[0] == "-format" {

		args = args[1:]

		max := len(args)
		if max < 1 {
			// no additional arguments, call original formatter
			ProcessFormat(doCompress, doCleanup)
			return
		}

		// first optional argument is parent pattern, will explore using Parent/* construct, write component XML with -element "*"
		prnt := args[0]
		if prnt == "" {
			fmt.Fprintf(os.Stderr, "\nPattern missing after -format command\n")
			os.Exit(1)
		}
		if prnt == "-xml" || prnt == "-doctype" || prnt == "-pfx" || prnt == "-sfx" {
			fmt.Fprintf(os.Stderr, "\nDeprecated argument '%s' used in -format command\n", prnt)
			os.Exit(1)
		}

		// second optional argument controls XML expansion or compression level
		elm := "*"
		addRet := false
		hideDoctype := false
		if max > 1 {
			// * = compact, ** = flush, *** = indented | @ = remove attributes | ^ = suppress xml and doctype
			numStars := 0
			hideAttrs := false
			for _, ch := range args[1] {
				if ch == '*' {
					numStars++
				} else if ch == '@' {
					hideAttrs = true
				} else if ch == '^' {
					hideDoctype = true
				}
			}
			if numStars > 1 {
				addRet = true
			}
			// construct legal element argument for PrintSubtree
			switch numStars {
			case 1:
				elm = "*"
			case 2:
				elm = "**"
			case 3:
				elm = "***"
			default:
				elm = "*"
			}
			if hideAttrs {
				elm += "@"
			}
		}

		fi, _ := os.Stdin.Stat()
		isPipe := bool((fi.Mode() & os.ModeCharDevice) == 0)

		if !isPipe {
			// no piped input, so write output instructions (without -head and -tail arguments)
			if addRet {
				fmt.Printf("xtract -pattern %s/* -ret \"\" -element \"%s\"\n", prnt, elm)
			} else {
				fmt.Printf("xtract -pattern %s/* -element \"%s\"\n", prnt, elm)
			}
			return
		}

		// add xml, DOCTYPE Parent, and <Parent> lines at the beginning
		hd := fmt.Sprintf("<?xml version=\"1.0\"?>\n<!DOCTYPE %s>\n<%s>", prnt, prnt)
		if hideDoctype {
			// or just <Parent> line
			hd = fmt.Sprintf("<%s>", prnt)
		}
		// add </Parent> line at the end
		tl := fmt.Sprintf("</%s>", prnt)

		// use -pattern Parent/* construct
		prnt += "/*"

		var acc []string

		acc = append(acc, "-head", hd, "-tail", tl)
		acc = append(acc, "-pattern", prnt)
		if addRet {
			acc = append(acc, "-ret", "")
		}
		acc = append(acc, "-element", elm)

		// data in pipe, so replace arguments, execute dynamically
		args = acc
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT

	head := ""
	tail := ""

	for {

		inSwitch = true

		switch args[0] {
		case "-head":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nPattern missing after -head command\n")
				os.Exit(1)
			}
			head = args[1]
		case "-tail":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nPattern missing after -tail command\n")
				os.Exit(1)
			}
			tail = args[1]
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past arguments
		args = args[2:]

		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "\nInsufficient command-line arguments supplied to xtract\n")
			os.Exit(1)
		}
	}

	// PARSE EXTRACTION ARGUMENTS

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nInsufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// check for -trial timing command
	trial := false
	if args[0] == "-trial" {
		trial = true
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nItem missing after -trial command\n")
			os.Exit(1)
		}
		if args[1] == "-pattern" || args[1] == "-Pattern" {
			// skip past -trial if followed by -pattern
			args = args[1:]
		} else {
			// otherwise substitute -pattern for -trial argument
			args[0] = "-pattern"
		}
	}

	// make sure top-level -pattern command is next
	if args[0] != "-pattern" && args[0] != "-Pattern" {
		fmt.Fprintf(os.Stderr, "\nNo -pattern in command-line arguments\n")
		os.Exit(1)
	}
	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nItem missing after -pattern command\n")
		os.Exit(1)
	}

	topPat := args[1]
	if topPat == "" {
		fmt.Fprintf(os.Stderr, "\nItem missing after -pattern command\n")
		os.Exit(1)
	}
	if strings.HasPrefix(topPat, "-") {
		fmt.Fprintf(os.Stderr, "\nMisplaced %s command\n", topPat)
		os.Exit(1)
	}

	// look for -pattern Parent/* construct for heterogeneous data, e.g., -pattern PubmedArticleSet/*
	topPattern, star := SplitInTwoAt(topPat, "/", LEFT)
	if topPattern == "" {
		return
	}

	parent := ""
	if star == "*" {
		parent = topPattern
	} else if star != "" {
		fmt.Fprintf(os.Stderr, "\n-pattern Parent/Child construct is not supported\n")
		os.Exit(1)
	}

	// parse nested exploration instruction from command-line arguments
	cmds := ParseArguments(args)
	if cmds == nil {
		fmt.Fprintf(os.Stderr, "\nProblem parsing command-line arguments\n")
		os.Exit(1)
	}

	// CREATE BUFFERED XML BLOCK READER

	in := bufio.NewReader(os.Stdin)
	rdr := NewXmlReader(in, doCompress, doCleanup)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create XML Block Reader\n")
		os.Exit(1)
	}

	// PERFORMANCE TIMING COMMAND

	// -trial tests speed of ProcessQuery streaming, requires valid extraction command
	if trial {

		// process directly, without using channels
		PartitionPattern(topPattern, star, rdr,
			func(rec int, str string) {
				ScanAndProcess(str, parent, cmds)
				recordCount++
				byteCount += len(str)
			})

		printDuration("records")

		return
	}

	// ORIGINAL PROCESSING METHOD IS RETAINED IN CASE OF PROBLEM WITH CONCURRENT SOLUTION

	// original logic where parser reads tokens directly from streamer
	if legacy {

		// group tokens by top pattern, parse and process each record
		recordCount = ProcessXml(cmds, topPattern, star, rdr)

		if timer {
			printDuration("records")
		}

		return
	}

	// LAUNCH PRODUCER AND CONSUMER SERVERS

	// chain processing steps into a concurrent pipeline with channels: producer -> consumers -> unshuffler -> results

	// launch producer goroutine to partition XML by pattern
	xmlq := CreateProducer(topPattern, star, chan_depth, rdr)
	if xmlq == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create producer\n")
		os.Exit(1)
	}

	// launch consumer goroutines to parse and explore partitioned XML objects
	tblq := CreateConsumers(cmds, parent, num_servers, chan_depth, xmlq)
	if tblq == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create consumers\n")
		os.Exit(1)
	}

	// launch unshuffler goroutine to restore order of results from autonomous consumers
	resq := CreateUnshuffler(chan_depth, heap_size, tblq)
	if resq == nil {
		fmt.Fprintf(os.Stderr, "\nUnable to create unshuffler\n")
		os.Exit(1)
	}

	// EXECUTE EXTRACTION COMMANDS

	first := true

	// drain reordered output channel
	for ext := range resq {
		str := ext.Text

		if str != "" {
			// on first non-empty result, print head, queue tail
			if first {
				first = false
				if head != "" {
					head = ConvertSlash(head)
					fmt.Printf("%s\n", head)
				}
				if tail != "" {
					tail = ConvertSlash(tail)
					defer fmt.Printf("%s\n", tail)
				}
			}
			if empty {
				// skip records with output
			} else if index {
				// print record index number
				idx := ext.Index
				fmt.Printf("%d\t%s", idx, str)
			} else {
				// normal output
				fmt.Printf("%s", str)
			}
		} else if empty {
			// report index numbers of records with no output
			idx := ext.Index
			fmt.Printf("%d\n", idx)
		}

		recordCount++

		runtime.Gosched()
	}

	if timer {
		printDuration("records")
	}
}
