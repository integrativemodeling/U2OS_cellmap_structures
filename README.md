# Structural models derived from the U2OS Cell Map 

## Summary

Human cells contain a complex hierarchy of components, many of which remain uncharted. Here we broadly identify subcellular components across the 10–8–10–5 meter range, via self-supervised fusion of immunofluorescence images and biophysical interactions measured for >5100 proteins in U-2 OS human cells. This global multiscale map resolves 275 molecular assemblies which we systematically validate by whole-cell size-exclusion chromatography and annotate using large language models. We explore key applications in driving structural biology, yielding structures for 111 heterodimeric complexes and an expanded 10-protein Rag-Ragulator assembly. The map assigns functions to 46 understudied proteins, including roles for C18orf21 in RNA processing and DPP9 in interferon signaling. It helps decode cancer genomes, identifying 21 recurrently mutated assemblies including cell-junction, nuclear pore, SWI/SNF and NCOR complexes, together implicating 102 cancer driver proteins that validate by transposon mutagenesis. The Multiscale Integrated Cell portal provides a queryable reference platform for structural and functional cell biology.


## List of files and directories:

### Rag Ragulator integrative structure model

- `data` All data used for integrative modeling, including the comparative model of the  SLC38A9-RagA-RagC-Ragulator complex, and AlphaFold2 predictions of BORCS6 and ITPA. 
- `scripts` PMI modeling script (`modeling.py`) to compute the integrative structure model of the Rag-Ragulator complex bound to BORCS6 and ITPA. 
- `analysis` Scripts to analyze the simulations 
- `results` All the relevant results from integrative modeling.
- `SI_table` Scripts to generate a table summarizing the integrative modeling protocols.
- `utils` Template and code to generate the <em>Supporting information</em> table summarizing the integrative modeling protocol.

### AlphaFold-Multimer Predictions

## Information

*Author (s)*: Ignacia Echeverria, Abantika Pal, Leah Schaffer, Clara Hu, Andrew Latham, Neelesh Soni 

_License_: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.

_Publications_: Multimodal cell maps as a foundation for structural and functional genomics. Leah V. Schaffer, Mengzhou Hu, Edward L. Huttlin, Gege Qian, Abantika Pal, Neelesh Soni, Andrew P. Latham, Laura Pontano Vaites, Kyung-Mee Moon, Dorothy Tsai, Nicole M. Mattson, Katherine Licon, Robin Bachelder, Anthony Cesnik , Ishan Gaur, Trang Le, William Leineweber, Aji Palar, Ernst Pulido , Yue Qin, Xiaoyu Zhao, Christopher Churas, Joanna Lenkiewicz, Jing Chen, Kei Ono, Dexter Pratt, Peter Zage, Ignacia Echeverria, Andrej Sali, Leonard J. Foster, J. Wade Harper, Steven P. Gygi, Emma Lundberg, Trey Ideker


