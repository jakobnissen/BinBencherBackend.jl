## [Input file formats](@id input_files)
### Binning TSV file
The binning is a tab-separated values (TSV) file of the following format:

* The header must be exactly `clustername\tcontigname`
* Every subsequent line must be the clustername, then a tab, then the contigname.

The format can be written as a regular expression:
```
binning = "clustername\tcontigname" ("\r?\n" bin)* ("\r?\n")?
bin = identifier "\t" identifier
identifier = "[^\t\r\n]*"
```

The file is assumed to be encoded in UTF-8. The identifiers may contain any arbitrary bytes
except '\t' (0x09), '\n' (0x0a) or '\r' (0x0d), including invalid unicode.

### [Reference JSON file](@id refjson)
This reference file is complex, as it contains a wealth of information about the genomes and their relationship.
It is recommended to create this file using the CLI - see [the relevant section of this documentation](@ref references).

If you want to create the JSON file yourself, you can learn its format by reading files in [the BinBencherBackend.jl repository](https://github.com/jakobnissen/BinBencherBackend.jl):
* A JSON Schema with descriptions at `files/schema.json`
* A small example reference in `files/ref.json`
* The source code in `src/reference.json` - grep for `ReferenceJSON`
