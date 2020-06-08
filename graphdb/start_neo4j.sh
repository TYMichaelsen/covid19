#!/bin/bash
# Current issues:  port map not working, ulimit warning persisting
cd /srv/rbd/covid19/bisystem
singularity instance start --fakeroot --writable-tmpfs --net --network-args "portmap=7474:7474/tcp" --bind neo4jdata:/var/lib/neo4j/data,neo4jrun:/var/lib/neo4j/run,neo4jlogs:/var/lib/neo4j/logs  neo4j.sif ne4j1
singularity exec --fakeroot instance://ne4j1 neo4j start
