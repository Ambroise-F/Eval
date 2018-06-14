#!/bin/bash

./evalverif > eval.txt
./ogverif > og.txt
diff eval.txt og.txt -s