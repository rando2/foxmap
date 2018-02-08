# Scripts for Chromosome Assembly of the Red Fox

These are the essential scripts developed in the assembly of the red fox genome. They are made available here in case they can be useful for future analyses. 

<b>essential_mods.py</b> is a set of very simple functions that we used to pull sequence from the scaffolds and orient it appropriately.<br> 
<b>racaout_to_csv.py</b> takes RACA's output files and transforms them to a csv format that is easier to use for analysis (in Excel, Python, R, etc.)
<b>make-frags.py</b> takes the scaffolds and the RACA input to make the fasta files that comprise 1) RACA's assembly, and 2) the chromosome fragment assembly. It depends on:<br>
<ol>
  <li>A comma-separated file (here: racafrag_foxfrag.csv) that specifies the relationship between RACA fragments and fox fragments. It shoudl have a header (Fragment,FoxSeg) and each line should list the RACA Fragment and its associated Chromosomal Fragment (e.g. 26a,10</li>
  <li>Syntenic relationships between the dog (reference) and fox (target) species, which are hard-coded into the script here as a list (line 113)</li>
  <li>the essential_mods.py functions</li>
  <li>RACA's output</li>
  <li>The scaffolds</li>
  </ol>

You can reach out to me at rando2 at Illinois dot edu if you have any questions.

MIT License

Copyright (c) [2018] [Halie M. Rando]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
