FLY = /data/Dfam/minifly

fly:
	rm -rf $(FLY)/sculu-out
	cargo run -- --logfile minifly.log \
		run \
		--alphabet   dna \
		--consensus  $(FLY)/consensi.fa \
		--alignments $(FLY)/seed_alignments \
		--config     $(FLY)/sculu.toml \
		--outdir     $(FLY)/sculu-out

run:
	cargo run -- --logfile - \
		run \
		--config                 sculu.toml \
		--alphabet               dna \
		--outdir                 out-test \
		--outfile                out-test/test.fa \
		--consensi               tests/inputs/consensi.fa \
		--instances              tests/inputs/instances 

build:
	cargo run -- --logfile - \
		components \
		--config                 sculu.toml \
		--alphabet               dna \
		--outdir                 out-test \
		--consensi               tests/inputs/consensi.fa \
		--instances              tests/inputs/instances

cluster:
	cargo run -- --logfile - \
		cluster \
		--config                 sculu.toml \
		--alphabet               dna \
		--outdir                 out-test \
		--consensi               out-test/consensi.fa \
		--instances              out-test/instances \
		--component              out-test/components/component-0

concat:
	cargo run -- --logfile - \
		concat \
		--outfile    out-test/families.fa \
		--consensi   out-test/consensi.fa \
		--components out-test/component-0/final.fa

mfs:
	cargo run -- --logfile   - \
		run \
		--config    sculu.toml \
		--alphabet  protein \
		--outdir    out-mfs \
		--consensi  ~/wheelerlab/data/proteins/mfs/consensi.fa \
		--instances ~/wheelerlab/data/proteins/mfs/instances 

pf:
	cargo run -- \
		--config    sculu.toml \
		--alphabet  protein \
		--outdir    pfam \
		--consensi  ~/wheelerlab/data/proteins/pfam/consensi.fa \
		--instances ~/wheelerlab/data/proteins/pfam/instances 

dup:
	cargo run -- \
		--config    sculu.toml \
		--instances tests/inputs/instances \
		--consensi  tests/inputs/dup_consensi.fa 

alu:
	cargo run -- --logfile - \
		run \
		--alphabet  dna \
		--config    sculu.toml \
		--consensi  ~/wheelerlab/data/alu/consensi.fa \
		--instances ~/wheelerlab/data/alu/instances \
		--outfile   ./alu-out/families.fa \
		--outdir    ./alu-out

mix:
	cargo run -- \
		--config    sculu.toml \
		--outdir    out-mixed \
		--consensi  data/mixed/consensi.fa \
		--instances data/mixed/instances \
		--outfile   final-mixed.fa \

tua:
	cargo run -- --logfile - \
		run \
		--alphabet  dna \
		--config    sculu.toml \
        --outdir    tuatara-out \
        --outfile   tuatara-out/new-consensi.fa \
        --instances ~/wheelerlab/tuatara/instances \
        --consensi  ~/wheelerlab/tuatara/consensi.fa

clean:
	rm -rf RM_* makedb.log
