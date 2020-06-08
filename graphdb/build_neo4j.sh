#!/usr/bin/env bash
cd ~
cp /srv/rbd/covid19/git/covid19/graphdb/neo4j.def .
singularity build --fakeroot neo4j.sif neo4j.def
cp neo4j.sif /srv/rbd/covid19/bisystem/
rm ~/neo4j.*
