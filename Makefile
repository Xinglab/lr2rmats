CC      =	gcc
CFLAGS  =	-Wall -O2 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  =	-g -Wall  
HTSLIB_DIR = ./htslib
HTSLIB  =   $(HTSLIB_DIR)/libhts.a
LIB     =	$(HTSLIB) -lm -lz -lpthread
INCLUDE = -I $(HTSLIB_DIR)

BIN_DIR =	./bin
SRC_DIR =   ./src

SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/lr2rmats
SORT 	=   sort_gtf.sh
GTF2GP  =   gtfToGenePred
GP2BED  =   genePredToBed

GDB_DEBUG   =   $(BIN_DIR)/gdb_lr2rmats
DMARCRO 	=	-D __DEBUG__

# dependencies
MINIMAP2  = minimap2
STAR      = STAR
SNAKEMAKE = snakemake
PSUTIL    = psutil
SAMTOOLS  = samtools
MINIMAP2_VERSION = 2.5
STAR_VERSION = 2.5.3a
SAMTOOLS_VERSION = 1.6

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all:		$(HTSLIB) $(BIN) $(GTF2GP) $(GP2BED)  $(SNAKEMAKE) $(SAMTOOLS) $(MINIMAP2) $(STAR)
lr2rmats_merge2names:     $(HTSLIB) $(BIN) $(GTF2GP) $(GP2BED)
lr2rmats:     $(HTSLIB) $(BIN) $(GTF2GP) $(GP2BED)
gdb_lr2rmats: $(SOURCE) $(GDB_DEBUG) 
dependencies: $(SNAKEMAKE) $(SAMTOOLS) $(MINIMAP2) $(STAR) 

$(SNAKEMAKE):
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	if [ -z ${shell which ${SNAKEMAKE}} ]; then \
		python3 -m pip install --user $(SNAKEMAKE) $(PSUTIL); \
		else echo "$(SNAKEMAKE) is already installed."; \
		fi

$(SAMTOOLS):
	if [ -z ${shell which ${SAMTOOLS}} ]; then \
		if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; \
		if [ ! -f ${BIN_DIR}/${SAMTOOLS} ]; then \
		wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2; \
		tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2; \
		cd ./samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}/ && ./configure --disable-bz2 --disable-lzma; make; \
		cd ..; ./configure --with-htslib=./htslib-${SAMTOOLS_VERSION}; make; \
		cp ${SAMTOOLS} ../${BIN_DIR}; cd .. ; \
		rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2 ./samtools-${SAMTOOLS_VERSION}; \
		else echo "$(SAMTOOLS) is already installed."; \
		fi; \
		else echo "$(SAMTOOLS) is already installed."; \
		fi

$(MINIMAP2):
	if [ -z ${shell which ${MINIMAP2}} ]; then \
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; \
		if [ ! -f ${BIN_DIR}/${MINIMAP2} ]; then \
		wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2; \
		tar -xjf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 || exit 255; \
		cp ./minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 ${BIN_DIR}; \
		rm -rf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 ./minimap2-${MINIMAP2_VERSION}_x64-linux; \
		else echo "$(MINIMAP2) is already installed."; \
		fi; \
		else echo "$(MINIMAP2) is already installed."; \
		fi

$(STAR):
	if [ -z ${shell which ${STAR}} ]; then \
		if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; \
		if [ ! -f ${BIN_DIR}/${STAR} ]; then \
		wget https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz; \
		tar -xzf ${STAR_VERSION}.tar.gz || exit 255; \
		cp ./STAR-${STAR_VERSION}/bin/Linux_x86_64/STAR ${BIN_DIR}; \
		rm -rf ${STAR_VERSION}.tar.gz ./STAR-${STAR_VERSION}; \
		else echo "$(STAR) is already installed."; \
		fi; \
		else echo "$(STAR) is already installed."; \
		fi

$(GTF2GP):
	if [ -z ${shell which ${GTF2GP}} ]; then \
		if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; \
		if [ ! -f ${BIN_DIR}/${GTF2GP} ]; then \
		wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/${GTF2GP}; \
		mv ./$(GTF2GP) ${BIN_DIR}; \
		else echo "$(GTF2GP) is already installed."; \
		fi; \
		else echo "$(GTF2GP) is already installed."; \
		fi

$(GP2BED):
	if [ -z ${shell which ${GP2BED}} ]; then \
		if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi; \
		if [ ! -f ${BIN_DIR}/${GP2BED} ]; then \
		wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/${GP2BED}; \
		mv ./$(GP2BED) ${BIN_DIR}; \
		else echo "$(GP2BED) is already installed."; \
		fi; \
		else echo "$(GP2BED) is already installed."; \
		fi


$(HTSLIB):
	cd $(HTSLIB_DIR); make;


$(BIN): $(OBJS)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(OBJS) -o $@ $(LIB)
	cp $(SRC_DIR)/$(SORT) $(BIN_DIR) 2> /dev/null


$(GDB_DEBUG):
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) $(INCLUDE) -o $@ $(LIB)

clean:
	rm -f $(SRC_DIR)/*.o $(BIN_DIR)/$(BIN) $(BIN_DIR)/$(SORT) 

clean_debug:
	rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(RGDB_DEBUG) $(NOR_DEBUG)
