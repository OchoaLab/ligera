# a fast way to remove correlated loci that are significant, keeping only most significant subset (ranked by p-value)
# simply returns indexes being kept
# - uses chr/pos info if available, not prune things that are far away, but just prune if it's missing
# TODO: document and export???
ld_prune <- function(
                     X,
                     pvals,
                     indexes,
                     r2_max = 0.3,
                     bim = NULL,
                     pos_window = 10000000L,
                     loci_on_cols = FALSE
                     ) {
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( pvals ) )
        stop( '`pvals` is required!' )
    if ( missing( indexes ) )
        stop( '`indexes` is required!' )

    # validate more (indexes should be actual numeric indexes, not booleans; bim should have pos and chr columns)
    

    # handle trivial singleton case
    if ( length( indexes ) == 1L )
        return( indexes )

    # override this for BEDMatrix
    if ( 'BEDMatrix' %in% class(X) ) {
        loci_on_cols <- TRUE
    } else if (!is.matrix(X))
        stop('X has unsupported class: ', toString( class( X ) ) )
    
    # order indexes by significance (original order unused/unimportant)
    indexes <- indexes[ order( pvals[ indexes ] ) ]

    # start navigating and deciding what to do
    # first one is always kept
    indexes_keep <- indexes[ 1L ]
    # remove from queue of things to test
    indexes <- indexes[ -1L ]
    # make sure top genotype isn't fixed (troubleshooting)
    # this should never happen in a real dataset
    x_test <- if ( loci_on_cols ) X[ , indexes_keep ] else X[ indexes_keep, ]
    if ( all( x_test == x_test[!is.na(x_test)][ 1L ], na.rm = TRUE ) )
        stop( 'Top locus was constant!' )
    while ( length( indexes ) > 1L ) {
        # separate an index to test
        index_test <- indexes[ 1L ]
        # and remove from queue
        indexes <- indexes[ -1L ]
        # get genotype vector now
        x_test <- if ( loci_on_cols ) X[ , index_test ] else X[ index_test, ]
        # an absurd edge case happens when one of the loci is constant!  In that case all its correlations are NA technically, but let's treat it as not significant to avoid straight out error, skipping entirely works out!
        # check if all values equal the first non-NA value
        if ( all( x_test == x_test[!is.na(x_test)][ 1L ], na.rm = TRUE ) )
            next
        # get coordinates if available
        if ( !is.null( bim ) ) {
            chr_test <- bim$chr[ index_test ]
            pos_test <- bim$pos[ index_test ]
        }
        # change if correlation fails
        keep_test <- TRUE
        
        # go in order and test correlation against every top locus
        for ( index_test2 in indexes_keep ) {
            # change this if bim is available and says otherwise
            check_pos <- TRUE
            if ( !is.null( bim ) && pos_window != 0L )
                # check optionally if chr/pos make sense for LD
                # has to be same chr and within window for check_pos to remain TRUE
                if ( chr_test != bim$chr[ index_test2 ] || abs( pos_test - bim$pos[ index_test2 ] ) > pos_window ) {
                    check_pos <- FALSE
            }

            # if !check_pos, then no correlation test is needed, these are treated as independent
            if ( check_pos ) {
                # actually calculate correlations here
                x_test2 <- if ( loci_on_cols ) X[ , index_test2 ] else X[ index_test2, ]
                # squared correlation!
                # this prevents errors if there are no complete cases, but it can cause warnings, which sadly we have to suppress for my tests to pass (and which I don't want my users to be troubled with).
                suppressWarnings(
                    r2_test <- stats::cor( x_test, x_test2, use = 'na.or.complete' )^2
                )
                # can be NA if entire missingness overlaps!
                # in that case, treat as independent! (only occurs in toy tests)
                if ( !is.na( r2_test ) && r2_test > r2_max ) {
                    # mark as not-keep
                    keep_test <- FALSE
                    # stop looking for more
                    break
                }
            }
        }

        # if not keeping, it was already removed from queue so nothing to do
        # otherwise add to list to keep, that's all that's needed
        if ( keep_test )
            indexes_keep <- c( indexes_keep, index_test )
    }

    # all done, return!
    return( indexes_keep )
}
