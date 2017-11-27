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

BIN     =	$(BIN_DIR)/lr2gtf
SORT 	=   sort_gtf.sh

GDB_DEBUG   =   $(BIN_DIR)/gdb_lr2gtf
DMARCRO 	=	-D __DEBUG__

# dependencies
MINIMAP2  = minimap2
STAR      = STAR
SNAKEMAKE = snakemake
SAMTOOLS  = samtools

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all:		$(HTS_LIB) $(BIN) 
lr2gtf:     $(BIN)
gdb_lr2gtf: $(SOURCE) $(GDB_DEBUG) 
dependendcies: $(SNAKEMAKE) $(SAMTOOLS) $(MINIMAP2) $(STAR) 

$(BIN): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LIB)
	cp $(SRC_DIR)/$(SORT) $(BIN_DIR) 2> /dev/null


$(SNAKEMAKE):
	if [ -z ${shell which ${SNAKEMAKE}} ]; then \
		git clone https://bitbucket.org/snakemake/snakemake.git && cd snakemake; \
		python3 setup.py install --user; cd ..; \
		rm -rf snakemake; \
		fi

$(SAMTOOLS):
	if [ ! -f ${BIN_DIR}/${SAMTOOLS} ]; then \
		wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2; \
		tar -jxvf samtools-1.6.tar.bz2; \
		cd ./samtools-1.6/ && make; cp ${SAMTOOLS} ../${BIN_DIR}; cd .. ; \
		rm -rf samtools-1.6.tar.bz2 ./samtools-1.6; \
		fi

$(MINIMAP2):
	if [ ! -f ${BIN_DIR}/${MINIMAP2} ]; then \
		wget https://github.com/lh3/minimap2/releases/download/v2.5/minimap2-2.5_x64-linux.tar.bz2; \
		tar -xjf minimap2-2.5_x64-linux.tar.bz2 || exit 255; \
		cp ./minimap2-2.5_x64-linux/minimap2 ${BIN_DIR}; \
		rm -rf minimap2-2.5_x64-linux.tar.bz2 ./minimap2-2.5_x64-linux; \
		fi

$(STAR):
	if [ ! -f ${BIN_DIR}/${STAR} ]; then \
		wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz; \
		tar -xzf 2.5.3a.tar.gz || exit 255; \
		cp ./STAR-2.5.3a/bin/Linux_x86_64/STAR ${BIN_DIR}; \
		rm -rf 2.5.3a.tar.gz ./STAR-2.5.3a; \
		fi


$(HTS_LIB):
	cd $(HTSLIB_DIR); make;

$(GDB_DEBUG):
	$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) $(INCLUDE) -o $@ $(LIB)

clean:
	rm -f $(SRC_DIR)/*.o $(BIN_DIR)/$(BIN) $(BIN_DIR)/$(SORT) 

clean_debug:
	rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(RGDB_DEBUG) $(NOR_DEBUG)
