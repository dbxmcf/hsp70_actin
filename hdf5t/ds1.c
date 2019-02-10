    /* Dataset1: each process takes a block of columns. */
    slab_set(start, count, stride, BYCOL);
if (verbose)
    printf("start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu\n",
    (unsigned long)start[0], (unsigned long)start[1],
        (unsigned long)count[0], (unsigned long)count[1],
        (unsigned long)(count[0]*count[1]));

    /* create a file dataspace independently */
    file_dataspace = H5Dget_space (dataset1);
    assert(file_dataspace != FAIL);
    MESG("H5Dget_space succeed");
    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride,
        count, NULL);
    assert(ret != FAIL);
    MESG("H5Sset_hyperslab succeed");

    /* create a memory dataspace independently */
    mem_dataspace = H5Screate_simple (SPACE1_RANK, count, NULL);
    assert (mem_dataspace != FAIL);

    /* fill dataset with test data */
    dataset_fill(start, count, stride, &data_origin1[0][0]);
    MESG("data_array initialized");
    if (verbose){
    MESG("data_array created");
    dataset_print(start, count, stride, &data_array1[0][0]);
    }

    /* set up the collective transfer properties list */
    xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist != FAIL);
    ret=H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(ret != FAIL);
    MESG("H5Pcreate xfer succeed");

    /* read data collectively */
    ret = H5Dread(dataset1, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
        xfer_plist, data_array1);
    assert(ret != FAIL);
    MESG("H5Dread succeed");

    /* verify the read data with original expected data */
    ret = dataset_vrfy(start, count, stride, &data_array1[0][0], &data_origin1[0][0]);
    assert(ret != FAIL);

    /* release all temporary handles. */
    /* Could have used them for dataset2 but it is cleaner */
    /* to create them again.*/
    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Pclose(xfer_plist);