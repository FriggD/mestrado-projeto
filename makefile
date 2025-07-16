passo-1:
	echo "Filtrando apenas os snps do arquivo do chr1"
	bcftools view -v snps ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o snps_chr1.vcf.gz

passo-2:
	echo "Criando o arquivo bed para depois calcular o ld"
	plink --vcf snps_chr1.vcf.gz --make-bed --out chr1

passo-3:
	echo "Calculando o ld" por janela
	plink --bfile chr1 --r2 --ld-window 99999 --ld-window-kb 500 --ld-window-r2 0.2 --out chr1_ld


download-genotypes:
	echo "Baixando de https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz.tbi


install-bcftools:
	sudo apt-get update
	sudo apt-get install gcc bzip2
	sudo apt-get install make
	sudo apt-get install libbz2-dev
	sudo apt-get install zlib1g-dev
	sudo apt-get install libncurses5-dev 
	sudo apt-get install libncursesw5-dev
	sudo apt-get install liblzma-dev

	echo "Instalando HTSLIB"
	cd /usr/bin
	sudo wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
	sudo tar -vxjf htslib-1.9.tar.bz2
	cd htslib-1.9
	sudo make

	echo "Instalado samtools"
	cd ..
	sudo wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	sudo tar -vxjf samtools-1.9.tar.bz2
	cd samtools-1.9
	sudo make

	echo "Instalando bcftools"
	cd ..
	sudo wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
	sudo tar -vxjf bcftools-1.9.tar.bz2
	cd bcftools-1.9
	sudo make

	export PATH="$PATH:/usr/bin/bcftools-1.9"
	export PATH="$PATH:/usr/bin/samtools-1.9"
	export PATH="$PATH:/usr/bin/htslib-1.9"
	source ~/.profile

install-plink:
	echo "https://knowledgebase.aridhia.io/workspaces/analysing-data/virtual-machines/installing-software-on-virtual-machines/installing-plink-on-your-virtual-machine"