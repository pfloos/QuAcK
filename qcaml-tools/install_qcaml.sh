#!/bin/bash

git clone --single-branch --branch dev https://gitlab.com/scemama/qcaml.git
cd qcaml
make
opam install .

