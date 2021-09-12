#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include "cgranges.h"
#include "htslib/sam.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

struct amp {
  char *name;
  int strand;
  int32_t st, en;
};

// parse a bed entry
char *parse_bed(char *s, int32_t *st_, int32_t *en_, struct amp *a_)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1, score;
	char *name;
	int strand;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) en = atol(p);
			else if (i == 3) name = strdup(p);
			else if (i == 4) score = atoi(p);
			else if (i == 5) strand = strcmp(p,"+")==0 ? 1 : 0;
			++i, p = q + 1;
			if (c == '\0') break;
		}
	}
	a_->name = name;
	a_->strand = strand;
	a_->st = st;
	a_->en = en;
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

// read bed file into a cgranges object and also store an array of structs w/ amplicon name and strand.

cgranges_t *read_bed(const char *fn, struct amp *a)
{
	gzFile fp;
	cgranges_t *cr;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	ks = ks_init(fp);
	cr = cr_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
	  char *ctg;
	  struct amp a1;
	  int32_t st, en;
	  ctg = parse_bed(str.s, &st, &en, &a1);
	  if (ctg) cr_add(cr, ctg, st, en, k);
	  a[k] = a1;
	  k++;
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return cr;
}

int main(int argc, char *argv[])
{

  int maxDist = 6;
  int verbose = 0;
  
 
  int c;
  opterr = 0;
  while ((c = getopt (argc, argv, "vd:")) != -1)
    switch (c)
      {
      case 'v':
        verbose = 1;
        break;
      case 'd':
        maxDist = atoi(optarg);
        break;
      case '?':
        if (optopt == 'd')
          fprintf (stderr, "Option -%c requires an integer argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
  
  if (argc - optind < 3) {
    printf("Usage: tagbam [options] <input.bam> <amplicon.bed> <output.bam>\n");
    printf("Options:\n");
    printf("  -d <int>   maximum sum of the distance bt amplicon and read ends (default=6)\n");
    printf("  -v         verbose\n");
    return 0;
  }

  char *inbam = argv[optind];
  char *inbed = argv[optind+1];
  char *outbam = argv[optind+2];
    
  // pointer array of structs to store amplicons
  struct amp *amplicons = (struct amp*)malloc(1000000*sizeof(struct amp));

  // read in bed file and store as a cgranges object + array of amplicons
  cgranges_t *cr;
  cr = read_bed(inbed,amplicons);
  assert(cr);
  cr_index(cr);

  fprintf(stderr,"loaded %d amplicon bed entries\n",(int)cr->n_r);  

  fprintf(stderr,"loading %s\n",inbam);
    
  // now open bam file
  samFile *fpin = sam_open(inbam, "rb");
  hts_idx_t *idx = sam_index_load(fpin, inbam);
  bam_hdr_t *h = sam_hdr_read(fpin);
  bam1_t *b = bam_init1();
  hts_itr_t *itr = sam_itr_querys(idx, h, ".");

  // make output bam file 
  samFile *fpout = sam_open(outbam, "wb");
  if (sam_hdr_write(fpout, h) != 0){
    fprintf(stderr,"cannot create bam file: %s",outbam);
    goto clean;
  }
  
  int readcounter = 0;
  int readstagged = 0;
  
  // iterate through bam file one at a time
  while (sam_itr_next(fpin, itr, b) >= 0){

    int strand;
    hts_pos_t st, en;
    const char *ref = h->target_name[b->core.tid];

    // dist = sum of distances between starts and ends of read and amplicon
    int dist = -1;
    // the index of hit that minimizes the sum of end distances
    int besthit = -1;

    // helper function to get end position of read
    if (!(b->core.flag & BAM_FREVERSE)){
      st = b->core.pos;
      en = st + llabs(b->core.isize) - 1;
      
    } else {
      st = b->core.mpos;
      en = st + llabs(b->core.isize) - 1;
    }
    
    // strand of read source molecule/probe: read1 is + = 1, read1 is - = 0 
    if (((b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREVERSE)) || ((b->core.flag & BAM_FREAD2) && (b->core.flag & BAM_FREVERSE))){
      strand = 1;
    } else {
      strand = 0;
    }      

    int numhits = 0;
    // only tag properly paired reads with amplicons
    if(b->core.flag & BAM_FPROPER_PAIR){
      
      int64_t i, n, *d = 0, max_d = 0;      
      n = cr_overlap(cr, ref, st, en, &d, &max_d);
      numhits = (int)n;
      for (i = 0; i < n; ++i){ // traverse overlapping intervals and find 
	if(amplicons[cr_label(cr, d[i])].strand != strand)
	  continue;

	if(amplicons[cr_label(cr, d[i])].st != cr_start(cr, d[i]) || amplicons[cr_label(cr, d[i])].en != cr_end(cr, d[i]))
	  fprintf(stderr,"coordinate mismatch...\n");
	
	// make this hit the best one if the end distance sum is smaller
	if(dist < 0 || dist > llabs(cr_start(cr, d[i]) - st) + llabs(cr_end(cr, d[i]) - en)){
	  dist = llabs(cr_start(cr, d[i]) - st) + llabs(cr_end(cr, d[i]) - en);
	  besthit = cr_label(cr, d[i]);
	}
      }
      // cleanup
      free(d);
           
      // if dist is <= maxDist then assign the read to the amplicon
      if (dist > 0 && dist <= maxDist){

	readstagged++;
	
	// delete the XN tag, if present
	uint8_t *data;
	if ((data = bam_aux_get(b,"XN")) != NULL)
	  bam_aux_del(b, data);

	bam_aux_append(b, "XN", 'Z', strlen(amplicons[besthit].name)+1, (uint8_t*)amplicons[besthit].name);
	  
      }
    }

    if(verbose==1){

      // this is for development purposes--read in the XN tag, if present, and then print out the read data.
      char *xn;
      uint8_t* dat;
      if ((dat = bam_aux_get(b,"XN")) != NULL)
	xn = bam_aux2Z(dat);

      fprintf(stderr,"%s\t%lld\t%lld\t%s\t%d\t%s\t%s\t%d\t%d\n", h->target_name[b->core.tid],st,en,bam_get_qname(b),strand,amplicons[besthit].name,xn,numhits,dist);
    }
    
    // write to new bam
    if (sam_write1(fpout, h, b) < 0){
      fprintf(stderr,"error writing to bam file\n");
      goto clean;
    }
  
    // progress counter
    readcounter++;
    if (readcounter % 10000 == 0)
      fprintf(stderr,"processed %d reads...\n",readcounter);    
  
  }

  fprintf(stderr,"tagged %d reads out of %d (%.2f%%) with maxDist=%d\n",readstagged,readcounter, (double)readstagged / (double)readcounter * 100.0, maxDist);
  
 clean:
  for (int i=0; i<1000000; i++) {
    free(amplicons[i].name);
  }
  free(amplicons);
  cr_destroy(cr);
  bam_destroy1(b);
  hts_itr_destroy(itr);
  bam_hdr_destroy(h);
  sam_close(fpin);
  sam_close(fpout);
  
}
