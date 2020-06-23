#!/usr/bin/env bash
cd ~
cp /srv/rbd/covid19/git/covid19/graphdb/neo4j_3_5.def .
singularity build --fakeroot neo4j_3_5.sif neo4j_3_5.def
cp neo4j_3_5.sif /srv/rbd/covid19/bisystem/
rm ~/neo4j.*
