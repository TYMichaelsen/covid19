#!/bin/bash
# Current issues:  port map not working, ulimit warning persisting
sudo singularity instance start --writable-tmpfs --net --network-args "portmap=7474:7474/tcp" --bind /data/neo4jdata:/var/lib/neo4j/data  neo4j.sif ne4j1
sudo singularity exec instance://ne4j1 neo4j start
