#!/bin/bash
wget 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL570&id=55999&db=GeoDb_blob143' -O GPL570-55999.txt
wget -i download_list.txt
