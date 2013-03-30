#! /usr/bin/env python
'''
Basic read/write tasks for mesh files.
'''
import numpy as np
import os
# ==============================================================================
def read(filenames, timestep=None):
    '''Reads an unstructured mesh with added data.

    :param filenames: The files to read from.
    :type filenames: str
    :param timestep: Time step to read from, in case of an Exodus input mesh.
    :type timestep: int, optional
    :returns mesh{2,3}d: The mesh data.
    :returns point_data: Point data read from file.
    :type point_data: dict
    :returns field_data: Field data read from file.
    :type field_data: dict
    '''
    if isinstance(filenames, (list, tuple)) and len(filenames)==1:
        filenames = filenames[0]

    if isinstance(filenames, basestring):
        filename = filenames
        # serial files
        extension = os.path.splitext(filename)[1]

        import re
        # setup the reader
        # TODO Most readers have CanReadFile() -- use that.
        if extension == '.vtu':
            from vtk import vtkXMLUnstructuredGridReader
            reader = vtkXMLUnstructuredGridReader()
            vtk_mesh = _read_vtk_mesh(reader, filename)
        elif extension == '.vtk':
            from vtk import vtkUnstructuredGridReader
            reader = vtkUnstructuredGridReader()
            vtk_mesh = _read_vtk_mesh(reader, filename)
        elif extension in [ '.ex2', '.exo', '.e' ]:
            from vtk import vtkExodusIIReader
            reader = vtkExodusIIReader()
            reader.SetFileName( filename )
            vtk_mesh = _read_exodusii_mesh(reader, timestep=timestep)
        elif re.match('[^\.]*\.e\.\d+\.\d+', filename):
            # Parallel Exodus files.
            # TODO handle with vtkPExodusIIReader
            from vtk import vtkExodusIIReader
            reader = vtkExodusIIReader()
            reader.SetFileName( filenames[0] )
            vtk_mesh = _read_exodusii_mesh(reader, timestep=timestep)
        else:
            raise RuntimeError( 'Unknown file type \'%s\'.' % filename )
    else:
        # Parallel files.
        # Assume Exodus format as we don't know anything else yet.
        from vtk import vtkPExodusIIReader
        # TODO Guess the file pattern or whatever.
        reader = vtkPExodusIIReader()
        reader.SetFileNames( filenames )
        vtk_mesh = _read_exodusii_mesh(reader, filename, timestep=timestep)

    return vtk_mesh
# ==============================================================================
def _read_vtk_mesh( reader, file_name ):
    '''Uses a vtkReader to return a vtkUnstructuredGrid.
    '''
    reader.SetFileName( file_name )
    reader.Update()

    return reader.GetOutput()
# ==============================================================================
def _read_exodusii_mesh( reader, timestep=None ):
    '''Uses a vtkExodusIIReader to return a vtkUnstructuredGrid.
    '''
    # Fetch metadata.
    reader.UpdateInformation()

    # Set time step to read.
    if timestep:
        reader.SetTimeStep( timestep )

    # Make sure the point fields are read during Update().
    for k in xrange( reader.GetNumberOfPointResultArrays() ):
        arr_name = reader.GetPointResultArrayName( k )
        reader.SetPointResultArrayStatus( arr_name, 1 )

    # Make sure all field data is read.
    for k in xrange( reader.GetNumberOfGlobalResultArrays() ):
        arr_name = reader.GetGlobalResultArrayName( k )
        reader.SetGlobalResultArrayStatus( arr_name, 1 )

    # Read the file.
    reader.Update()
    out = reader.GetOutput()

    # Loop through the blocks and search for a vtkUnstructuredGrid.
    vtk_mesh = []
    for i in xrange( out.GetNumberOfBlocks() ):
        blk = out.GetBlock( i )
        for j in xrange( blk.GetNumberOfBlocks() ):
            sub_block = blk.GetBlock( j )
            if sub_block.IsA( 'vtkUnstructuredGrid' ):
                vtk_mesh.append( sub_block )

    if len(vtk_mesh) == 0:
        raise IOError( 'No \'vtkUnstructuredGrid\' found!' )
    elif len(vtk_mesh) > 1:
        raise IOError( 'More than one \'vtkUnstructuredGrid\' found!' )

    # Cut off trailing '_' from array names.
    for k in xrange( vtk_mesh[0].GetPointData().GetNumberOfArrays() ):
        array = vtk_mesh[0].GetPointData().GetArray(k)
        array_name = array.GetName()
        if array_name[-1] == '_':
            array.SetName( array_name[0:-1] )

    return vtk_mesh[0]
# ==============================================================================
def write(
          filename,
          vtk_mesh,
          point_data = None,
          cell_data = None,
          field_data = None
          ):
    '''Writes mesh together with data to a file.

    :params filename: File to write to.
    :type filename: str

    :params point_data: Named additional point data to write to the file.
    :type point_data: dict
    '''
    import os

    extension = os.path.splitext(filename)[1]
    # add point data
    is_exodus_format = extension in [ '.ex2', '.exo', '.e' ]
    if point_data:
        for name, data in point_data.iteritems():
            new_name = name
            # There is a naming inconsistency in VTK when it comes to
            # multivectors in Exodus files:
            # If a vector 'v' has two components, they are called 'v_r',
            # 'v_z' (note the underscore), if it has three, then they are
            # called 'vx', 'vy', 'vz'.
            # Make this consistent by appending an underscore if needed.
            # Note that for VTK files, this problem does not occur since
            # the label of a vector is always stored as a string.
            is_3d_vector = len(data.shape) == 2 and data.shape[1] == 3
            if is_exodus_format and is_3d_vector and name[-1] != '_':
                new_name += '_'
            vtk_mesh.GetPointData() \
                    .AddArray(_create_vtkarray(data, new_name))

    # add cell data
    if cell_data:
        for key, value in cell_data.iteritems():
            vtk_mesh.GetCellData() \
                    .AddArray(_create_vtkarray(value, key))

    # add field data
    if field_data:
        for key, value in field_data.iteritems():
            vtk_mesh.GetFieldData() \
                    .AddArray(_create_vtkarray(value, key))

    import re
    extension = os.path.splitext(filename)[1]
    if extension == '.vtu': # VTK XML format
        from vtk import vtkXMLUnstructuredGridWriter
        writer = vtkXMLUnstructuredGridWriter()
    elif extension == '.pvtu': # parallel VTK XML format
        from vtk import vtkXMLPUnstructuredGridWriter
        writer = vtkXMLPUnstructuredGridWriter()
    elif extension == '.vtk': # classical VTK format
        from vtk import vtkUnstructuredGridWriter
        writer = vtkUnstructuredGridWriter()
        writer.SetFileTypeToASCII()
    elif extension in [ '.ex2', '.exo', '.e' ]: # Exodus II format
        from vtk import vtkExodusIIWriter
        writer = vtkExodusIIWriter()
        # If the mesh contains vtkModelData information, make use of it
        # and write out all time steps.
        writer.WriteAllTimeStepsOn()
    elif re.match('[^\.]*\.e\.\d+\.\d+', filename):
        # TODO handle parallel I/O with vtkPExodusIIWriter
        from vtk import vtkExodusIIWriter
        writer = vtkExodusIIWriter()
        # If the mesh contains vtkModelData information, make use of it
        # and write out all time steps.
        writer.WriteAllTimeStepsOn()
    else:
        raise IOError( 'Unknown file type \'%s\'.' % filename )

    writer.SetFileName( filename )

    writer.SetInput( vtk_mesh )

    writer.Write()

    return
# ==============================================================================
def _create_vtkarray(X, name):
    import numpy as np
    from vtk import vtkBitArray, vtkIntArray, vtkDoubleArray, vtkCharArray

    # If something isn't a Numpy array already, try to make it one.
    if not isinstance(X, np.ndarray) and not isinstance(X, str):
        X = np.array(X)

    # This could be a lot more fine-grained:
    # vtkLongLongArray, vtkFloatArray,...
    if isinstance(X, str) or X.dtype.kind == 'S':
        array = vtkCharArray()
    elif X.dtype == bool:
        array = vtkBitArray()
    elif X.dtype == int:
        array = vtkIntArray()
    elif X.dtype == float:
        array = vtkDoubleArray()
    elif X.dtype == complex:
        # Convert complex arrays to double.
        Y = np.empty((len(X),2), dtype=float)
        if len(X.shape) == 1:
            Y[:,0] = X.real
            Y[:,1] = X.imag
        elif len(X.shape) == 2:
            Y[:,0] = X[:,0].real
            Y[:,1] = X[:,0].imag
        else:
            raise RuntimeError()
        X = Y
        array = vtkDoubleArray()
    else:
        raise TypeError('Unknown VTK data type', X.dtype, '.')

    # For some reason, setting the number of tuples and then using
    # SetNextTuple() or similar doesn't seem to work:
    # The application segfaults or, worse, yields an irrecoverable
    # glibc memory corruption.
    # Most likely the cause: You have to call SetNumberOfTuples()
    # AFTER SetNumberOfComponents().
    #array.SetNumberOfTuples(X.shape[0])
    # Special treatment for strings:
    if isinstance(X, str):
        array.SetNumberOfComponents(len(X))
        array.SetNumberOfTuples(1)
        array.SetTupleValue(0, X)
    elif len(X.shape) == 0:
        array.SetNumberOfComponents(1)
        # Set values.
        array.InsertNextValue(X)
    elif len(X.shape) == 1:
        array.SetNumberOfComponents(1)
        # Set values.
        for k in xrange(X.shape[0]):
            array.InsertNextValue(X[k])
    elif len(X.shape) == 2:
        array.SetNumberOfComponents(X.shape[1])
        # Set values.
        for k in xrange(X.shape[0]):
            for k2 in xrange(X.shape[1]):
                array.InsertNextValue(X[k][k2])
    else:
        raise ValueError('Don''t know what to do with many-dimensional array ''%s''.' % name)

    array.SetName(name)

    return array
# ==============================================================================
