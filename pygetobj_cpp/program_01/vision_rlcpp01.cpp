#include "Python.h"
#include "c24bitmap.h"
#include <numpy/arrayobject.h>

/* Call double avg(double *, int) */
static PyObject *py_imgcpt(PyObject *self, PyObject *args) {
  PyObject *bufobj;
  Py_buffer view;
  double result;
  /* Get the passed Python object */
  if (!PyArg_ParseTuple(args, "O", &bufobj)) {
    return NULL;
  }

  /* Attempt to extract buffer information from it */
  if (PyObject_GetBuffer(bufobj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1) {
    return NULL;
  }

  printf("========================\n",view.ndim);
  printf("the opencv obj's dimensions are %i\n",view.ndim);
  printf("========================\n",view.ndim);

  printf("============ok0=============\n");
  printf(" shape:== %d, %d, %d =\n",    view.shape[0], view.shape[1], view.shape[2]);
  printf("stride:== %d, %d, %d =\n",  view.strides[0], view.strides[1], view.strides[2]);
  printf("====data format:%s=====\n",view.format);
  printf("============ok1=============\n");
  
  C24BitMap CPic;
  CPic.FormatF(view.shape[1], view.shape[0]);
  BYTE* ImgPtr = (BYTE*)view.buf;
 
  int i,j;

  Loopi(view.shape[0])
    { 
      BYTE * LineStart = CPic.GetRowStartPt(i);
      BYTE * DataSrc   = &ImgPtr[i*view.strides[0]];
     Loopj(view.shape[1]*3)
     {
       LineStart[j] =  DataSrc[j];
     }/* */
   }

   CPic.FillRect = 0;
   CPic.SetPenColor(255,0,0);
   CPic.Rectangle(10,10,400,300);

   CPic.SetPenColor(0,0,255);
   CPic.Rectangle(20,20,350,250);

   CPic.SetPenColor(0,255,0);
   CPic.Rectangle(30,30,300,200);

   CPic.Save("draw_temp.bmp");


  /* Pass the raw buffer and size to the C function */
  result = 10;//imgcpt(view.buf, view.shape[0]);

  /* Indicate we're done working with the buffer */
  PyBuffer_Release(&view);
  return Py_BuildValue("d", result);
}

/* Module method table */
static PyMethodDef SampleMethods[] = {
  {"imgcpt",  py_imgcpt, METH_VARARGS, "cpp_ext_calculation"},
  { NULL, NULL, 0, NULL}
};

/* Module structure */
static struct PyModuleDef samplemodule = {
  PyModuleDef_HEAD_INIT,
  "sample1",           /* name of module */
  "A sample module",  /* Doc string (may be NULL) */
  -1,                 /* Size of per-interpreter state or -1 */
  SampleMethods       /* Method table */
};

/* Module initialization function */
PyMODINIT_FUNC
PyInit_sample1(void) {
  return PyModule_Create(&samplemodule);
}
