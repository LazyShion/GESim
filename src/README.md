## GESim
The implementation of an entropy-based graph similarity search for a graph database.

## How to Use
### Requirements
This software requires the following softwares.
* gcc Version 11.3.1 (or later)

We have confirmed that our software works on the following environments.
* CentOS Stream release 9
* MacOS 12.6

### Build
After installing gcc, you can build this software as follows:
```
$ make clean; make
```
If you find `search` and `convert`, the build has been successfully completed.


### How to generate an input file from SMILES
#### Input SMILES file
Prepare an input SMILES file like as follows:

```
0 CCN(c1ccccc1)c1cc(NC)[n+](C)c(C)n1
1 CCOc1ccccc1C(=O)NCC1(N2CCN(CC)CC2)CCCCC1
2 CCOc1ccccc1C(=O)NCC1(N2CCN(C3CC3)CC2)CCCCC1
3 CCOc1ccccc1C(=O)NCC1(N2CCN(C(C)C)CC2)CCCCC1
4 COc1ccccc1C(=O)NCC1(N2CCN(C(C)C)CC2)CCCCC1
5 CCOc1ccccc1C(=O)NCC1(N2CCN(C3CCCC3)CC2)CCCCC1
6 CCOc1ccccc1C(=O)NCC1(N2CCN(C3CCC3)CC2)CCCCC1
7 CCOc1ccccc1C(=O)NCC1(N2CCN(C)CC2)CCCCC1
```

#### Generate a graph file from SMILES file
Run the following python script to generate a graph file from a SMILES file.
```
$ python smiles2graph.py <SMILES file> <output file name> <mapping file name>
```


### Usage


#### Input file
Input file need to include all graphs in a graph. The input file is specifically formatted as follows:

``` sample_graphDB.txt
t 0
v 0 C
v 1 O
v 2 C
...
v 38 C
v 39 C
e 0 0 1 1.0
e 1 1 2 1.0
e 2 2 3 1.5
e 3 3 4 1.5
e 4 4 5 1.0
e 5 5 6 1.0
e 6 4 7 1.5
e 7 7 8 1.5
...
e 41 30 9 1.5
e 42 39 2 1.5
t 1
v 0 C
v 1 O
v 2 C
...
```
`sample_graphDB.txt` includes 1,000 weighted labeld graphs.
One graph is descripted in the order of a graph-ID `t`, node-IDs `v`, and edge-IDs `e`, which are detailed as follows: 

|ID           |Syntax        |Semantics                               |
|-------------|:-------------|:---------------------------------------|
|Graph-ID `t` | `t <TAB> (int)n` |This line assigns node-ID `n` to a graph|
|Node-ID `v`  | `v <TAB> (int)n <TAB> (String)l` | This line descripts a node whose node-ID is `n` having a node-label `l`.|
|Edge-ID `e`  | `e <TAB> (int)n <TAB> (int)s <TAB> (int)d <TAB> (float)w` | This line descripts an edge of edge-ID `n` connecting node-IDs `s` and `d` with an edge weight `w`|

#### File conversion
This software reads a given graph DB by the CRS format, and this requires a file conversion process. 
To covert the input file (`sample_graphDB.txt`) into the CRS format (`sample_graphDB.bin`), you should run `convert` like as follows:
``` convert
$ ./convert -i sample_graphDB.txt -o sample_graphDB.bin
```
`convert` requires two options, `-i` and `-o`, that specify names of the input file and the CRS formatted file, respectively.

#### Query processing
Finally, we can run queries on the database by using `search` like as follows:
```
$ ./search -i sample_graphDB.bin -q 0
```
`search` has the following two options.
|Option|Description|
|------|:----------|
|`-i`  |Specify an input graph database (i.e., `sample_graphDB.bin`). |
|`-q`  |Specify a query node-ID in the database. |
|`-d` (optional)  |Specify a maximum diameter to build fingerprints. |

