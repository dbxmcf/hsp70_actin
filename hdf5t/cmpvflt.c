/*
   cmpvflt.c
 */

#include "hdf5.h"
#include "malloc.h"

typedef struct s1 {
   int type;
   int nElements;
   hvl_t pElements[1];
   float max;
} s1;
s1  st[3];

int main(int argc, char* argv[])
{

    int i,j, len, nCount  = 2;
    hsize_t     dims[1], chunkDims[1], ioffset[1], icnt[1], arr_dims[]={1};
    herr_t      status, err1, err2;
    hid_t fileId, dataSetId, dataSpaceId, memSpaceId, vlDataTypeId, dataTypeId, pListId;
    hid_t  array_dt, vltid;


    /* Create the file */
    fileId = H5Fcreate("cmpvflt.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    printf ("H5Fcreate returns: %i\n", fileId);

    /* Create the array datatype */


    vlDataTypeId= H5Tvlen_create ( H5T_NATIVE_FLOAT);

    array_dt = H5Tarray_create (vlDataTypeId, 1, arr_dims, NULL); 
    dataTypeId = H5Tcreate(H5T_COMPOUND, sizeof(s1));
    printf ("H5Tcreate returns: %i\n", dataTypeId);
    status = H5Tinsert(dataTypeId, "type", HOFFSET(s1,type), H5T_NATIVE_INT);
    status = H5Tinsert(dataTypeId, "nElements", HOFFSET(s1,nElements), H5T_NATIVE_INT);
    status = H5Tinsert(dataTypeId, "pElements", HOFFSET(s1,pElements), array_dt);
    status = H5Tinsert(dataTypeId, "max", HOFFSET(s1,max), H5T_NATIVE_FLOAT);

    /* Create the dataspace */
    dims[0] = nCount;
    dataSpaceId = H5Screate_simple(1, dims, NULL);
    printf ("H5Screate_simple returns: %i\n", dataSpaceId);

    /* Create the property list  */
    pListId = H5Pcreate(H5P_DATASET_CREATE);

    /* Create the dataset */
    chunkDims[0] = 1;
    status = H5Pset_chunk(pListId, 1, chunkDims);
    printf ("%i\n", status);
    dataSetId = H5Dcreate(fileId, "/SAMPLE4", dataTypeId, dataSpaceId, pListId);

    for (i=0; i<nCount; i++)
    {
	st[i].type = 1;
        len = (i+1)*10;
        st[i].nElements = len;
        st[i].pElements[0].len = len;
        st[i].pElements[0].p = malloc (len* sizeof (float));
        for (j=0; j<len; j++)
              ((float *)st[i].pElements[0].p)[j] = .5+j;
        st[i].max = 3.1;
    }

    err2 = H5Dwrite (dataSetId, dataTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, &st);
    printf ("H5Dwrite returns: %i\n", err2);

    status = H5Dvlen_reclaim(dataTypeId, dataSpaceId, H5P_DEFAULT, &st);

    /* Close everything */
    status= H5Pclose(pListId);
    printf ("%i\n", status);
    status = H5Dclose(dataSetId);
    printf ("status: %i\n", status);
    status = H5Sclose(dataSpaceId);
    printf ("status: %i\n", status);
    status = H5Tclose(vlDataTypeId);
    printf ("status: %i\n", status);
    status = H5Tclose(array_dt);
    printf ("status: %i\n", status);
    status = H5Tclose(dataTypeId);
    printf ("status: %i\n", status);
    status = H5Fclose(fileId);    
    printf ("status: %i\n", status);

    return 0;
}



