#include <Rcpp.h>
#include <cerrno>
#include <vector>
using namespace Rcpp;

// expected header (magic numbers)
// assume standard locus-major order and latest format
const unsigned char plink_bed_byte_header[3] = {0x6c, 0x1b, 1};

// this function calculates two terms from one pass of the genotype matrix, needed for unbiased implicit kinship estimates

// [[Rcpp::export]]
List get_b_inbr_bed_cpp(
			const char* file,
			std::vector<int>::size_type m_loci,
			std::vector<int>::size_type n_ind,
			double mean_kinship
			) {
  // file must be full path (no missing extensions)
  
  // open input file in "binary" mode
  FILE *file_stream = fopen( file, "rb" );
  // die right away if needed, before initializing buffers etc
  if ( file_stream == NULL ) {
    // send error message to R
    char msg[100];
    sprintf(msg, "Could not open BED file `%s` for reading: %s", file, strerror( errno ));
    stop(msg);
  }

  /////////////////////////
  // check magic numbers //
  /////////////////////////

  // for header only
  unsigned char buffer_header[3];
  // for extra sanity checks, keep track of bytes actually read (to recognize truncated files)
  // reuse this one for genotypes below
  std::vector<int>::size_type n_buf_read;
  
  // read header bytes (magic numbers)
  n_buf_read = fread( buffer_header, sizeof(unsigned char), 3, file_stream );
  // this might just indicate an empty file
  if ( n_buf_read != 3 ) {
    // wrap up everything properly
    fclose( file_stream ); // close file
    // now send error message to R
    stop("Input BED file did not have a complete header (3-byte magic numbers)!");
  }
  
  // require that they match our only supported specification of locus-major order and latest format
  // was using strcmp but there are funky issues (wants signed, but we don't really want order anyway, just test for equality)
  // use explicit loop instead
  int pos;
  for (pos = 0; pos < 3; pos++) {
    if ( plink_bed_byte_header[pos] != buffer_header[pos] ) {
      // wrap up everything properly
      fclose( file_stream ); // close file
      // now send error message to R
      stop("Input BED file is not in supported format.  Either magic numbers do not match, or requested sample-major format is not supported.  Only latest locus-major format is supported!");
    }
  }

  // preallocate the things we want:
  
  // (1) a summary of the MAF variance across the genome
  double b = 0.0;
  // intermediates for b, to calculate allele frequencies
  int x_sum = 0;
  std::vector<int>::size_type x_num = n_ind; // decrement NAs, just as for inbr_num above (assuming NAs are rare, this is faster)
  double x_mean;
  
  // (2) intermeditates for inbr
  // an uncorrected inbreeding coefficient estimate
  std::vector<int>::size_type* inbr_sum = new std::vector<int>::size_type[ n_ind ];
  // the number of non-NA cases, per individual
  std::vector<int>::size_type* inbr_num = new std::vector<int>::size_type[ n_ind ];
  std::vector<int>::size_type j;
  // initialize
  for (j = 0; j < n_ind; j++) {
    inbr_sum[ j ] = 0;
    // all start as m_loci, will decrement later as we encounter NAs
    inbr_num[ j ] = m_loci;
  }
  
  ////////////////////
  // read genotypes //
  ////////////////////
  
  // number of columns (bytes) in input (for buffer), after byte compression
  // size set for full row, but overloaded used first for this header comparison
  // chose std::vector<int>::size_type to have it match n_buf_read value returned by fread
  std::vector<int>::size_type n_buf = ( n_ind + 3 ) / 4;
  // initialize row buffer
  unsigned char* buffer = new unsigned char[ n_buf ];
  
  // navigate data and process
  std::vector<int>::size_type i;
  std::vector<int>::size_type k; // to match n_buf type
  unsigned char buf_k; // working of buffer at k'th position
  unsigned char xij; // copy of extracted genotype
  for (i = 0; i < m_loci; i++) {
    
    // read whole row into buffer
    n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
    
    // always check that file was not done too early
    if ( n_buf_read != n_buf ) {
      // wrap up everything properly
      delete [] buffer; // free buffer memory
      fclose( file_stream ); // close file
      // now send error message to R
      char msg[100];
      sprintf(msg, "Truncated file: row %ld terminated at %ld bytes, expected %ld.", i+1, n_buf_read, n_buf); // convert to 1-based coordinates
      stop(msg);
    }

    // process buffer now!

    // always reset these at start of row
    j = 0; // individuals

    // reset things we need to get allele frequency
    x_sum = 0;
    x_num = n_ind; // decrement NAs, just as for inbr_num above (assuming NAs are rare, this is faster)
      
    // navigate buffer positions k (not individuals j)
    for (k = 0; k < n_buf; k++) {
      
      // copy down this value, which will be getting edited
      buf_k = buffer[k];

      // navigate the four positions
      // pos is just a dummy counter not really used except to know when to stop
      // update j too, accordingly
      for (pos = 0; pos < 4; pos++, j++) {

	if (j < n_ind) {
	  // extract current genotype using this mask
	  // (3 == 00000011 in binary)
	  xij = buf_k & 3;
	
	  // handle cases (in original encoding; canonical encoding in comments)
	  // Homozygotes are more common, so test for those first
	  // inbr_sum is incremented upon each homozygote only!
	  // x_sum increments by ordinary xij, no need to overwrite xij at all
	  if (xij == 0) {
	    // xij = 2; // 0 -> 2
	    x_sum += 2;
	    inbr_sum[ j ]++;
	  } else if (xij == 3) {
	    // xij = 0; // 3 -> 0
	    inbr_sum[ j ]++;
	  } else if (xij == 2) {
	    // xij = 1; // 2 -> 1
	    x_sum += 1;
	  } else { // only case left is NA
	    // xij = 0; // 1 -> NA
	    // decrement non-NA count for this individual and locus, respectively
	    inbr_num[ j ]--;
	    x_num--;
	  }
	  
	  // shift packed data, throwing away genotype we just processed
	  buf_k = buf_k >> 2;
	} else {
	  // when j is out of range, we're in the padding data now
	  // as an extra sanity check, the remaining data should be all zero (that's how the encoding is supposed to work)
	  // non-zero values would strongly suggest that n_ind was not set correctly
	  if (buf_k != 0) {
	    // wrap up everything properly
	    delete [] buffer; // free buffer memory
	    fclose( file_stream ); // close file
	    // now send error message to R
	    char msg[200];
	    sprintf(msg, "Row %ld padding was non-zero.  Either the specified number of individuals is incorrect or the input file is corrupt!", i+1); // convert to 1-based coordinates
	    stop(msg);
	  }
	}
      }
      // finished byte
      
    }
    // finished row

    // finish processing terms for b
    // calculate mean genotype now that we're done with this one locus
    x_mean = static_cast<double>(x_sum) / x_num;
    // update b running sum
    b += x_mean * (2 - x_mean);
    
  }
  // finished matrix!

  // let's check that file was indeed done
  n_buf_read = fread( buffer, sizeof(unsigned char), n_buf, file_stream );
  // wrap up regardless
  delete [] buffer; // free buffer memory
  // and more troubleshooting messages (for windows)
  if ( fclose( file_stream ) != 0 )
    stop("Input BED file stream close failed!");
  if ( n_buf_read != 0 ) {
    // now send error message to R
    stop("Input BED file continued after all requested rows were read!  Either the specified the number of loci was too low or the input file is corrupt!");
  }

  // do some final processing of b
  // this includes all final steps, including the mean_kinship correction
  b = ( 1.0 - ( b / m_loci ) - mean_kinship ) / ( 1.0 - mean_kinship );
  // create R version
  NumericVector bR(1);
  bR[ 0 ] = b;

  // final inbreeding vector
  NumericVector inbr(n_ind);
  // the inbreeding data also now gets normalized fully, which requires the b above
  for (j = 0; j < n_ind; j++) {
    inbr[ j ] = ( ( 2.0 * inbr_sum[ j ] ) / inbr_num[ j ] - 1.0 - b ) / ( 1.0 - b );
  }

  // free intermediate objects now
  delete [] inbr_sum;
  delete [] inbr_num;

  // should return b and inbr together in an R List
  List ret_list = List::create(
			       Named("b") = bR,
			       Named("inbr") = inbr
			       );
  
  // return genotype matrix
  return ret_list;
}
