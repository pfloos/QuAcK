#include <Python.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>


void call_mol_prop(const char *xyz, const char *input_basis, int charge, int multiplicity, 
                   const char *unit, int cartesian, 
                   int *natoms, int *nalpha, int *nbeta, double *Enuc) {

    PyObject *pName, *pModule, *pFunc, *pArgs, *pValue;
    PyObject *py_xyz, *py_input_basis, *py_unit, *py_cartesian;
    PyObject *py_result;
    
    // Initialize the Python interpreter
    Py_Initialize();

    //printf("Python version: %s\n", Py_GetVersion());
    //printf("xyz: %s\n", xyz);
    //printf("input_basis: %s\n", input_basis);
    //printf("charge: %d\n", charge);
    //printf("multiplicity: %d\n", multiplicity);
    //printf("unit: %s\n", unit);
    //printf("cartesian: %d\n", cartesian);

    char exe_path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
    if (len != -1) {
        exe_path[len] = '\0'; // Null-terminate the string
        dirname(exe_path); // Get directory part
        char sys_path_command[PATH_MAX + 30];
        snprintf(sys_path_command, sizeof(sys_path_command), "sys.path.append('%s')", exe_path);
        PyRun_SimpleString("import sys");
        PyRun_SimpleString(sys_path_command);
    } else {
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append('.')");
    }
    pName = PyUnicode_DecodeFSDefault("pyscf_module");
    pModule = PyImport_Import(pName);
    Py_XDECREF(pName);

    if (pModule != NULL) {

        // Get the Python function
        pFunc = PyObject_GetAttrString(pModule, "mol_prop");

        if (pFunc && PyCallable_Check(pFunc)) {

            // Convert C strings to Python strings
            py_xyz = PyUnicode_FromString(xyz);
            py_input_basis = PyUnicode_FromString(input_basis);
            py_unit = PyUnicode_FromString(unit);

            // Convert C int to Python boolean
            py_cartesian = PyBool_FromLong(cartesian);

            // Create a tuple to hold the arguments for the Python function
            pArgs = PyTuple_Pack(6, 
                                 py_xyz, 
                                 py_input_basis, 
                                 PyLong_FromLong(charge), 
                                 PyLong_FromLong(multiplicity), 
                                 py_unit, 
                                 py_cartesian);
            
            // Call the Python function
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_XDECREF(pArgs);

            if (pValue != NULL) {
                // Unpack the result tuple
                if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 4) {
                    *natoms = (int)PyLong_AsLong(PyTuple_GetItem(pValue, 0));
                    *nalpha = (int)PyLong_AsLong(PyTuple_GetItem(pValue, 1));
                    *nbeta = (int)PyLong_AsLong(PyTuple_GetItem(pValue, 2));
                    *Enuc = PyFloat_AsDouble(PyTuple_GetItem(pValue, 3));
                } else {
                    PyErr_Print();
                }
                Py_XDECREF(pValue);  // Clean up reference to result tuple
            } else {
                // Handle error if Python function call fails
                PyErr_Print();
            }

            // Clean up references
            Py_XDECREF(py_xyz);
            Py_XDECREF(py_input_basis);
            Py_XDECREF(py_unit);
            Py_XDECREF(py_cartesian);
        } else {
            // Handle error if Python function is not callable
            PyErr_Print();
        }

        Py_XDECREF(pFunc);  // Clean up reference to Python function object
        Py_XDECREF(pModule);  // Clean up reference to Python module object
    } else {
        // Handle error if Python module could not be imported
        PyErr_Print();
    }

    // Finalize the Python interpreter
    Py_Finalize();
}

