////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "vtkWriter.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedIntArray.h>

rw::writer::VtkWriter::VtkWriter(const std::string &filename,
                                 const std::string &compress_type)
    : d_compressType(compress_type) {
  std::string f = filename + ".vtu";

  d_writer_p = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  d_writer_p->SetFileName(const_cast<char *>(f.c_str()));
}

void rw::writer::VtkWriter::appendNodes(const std::vector<util::Point3> *nodes,
                                        const std::vector<util::Point3> *u) {
  auto points = vtkSmartPointer<vtkPoints>::New();

  for (size_t i = 0; i < nodes->size(); i++) {
    util::Point3 p = (*nodes)[i];
    if (u) p = p + (*u)[i];
    points->InsertNextPoint(p.d_x, p.d_y, p.d_z);
  }

  d_grid_p = vtkSmartPointer<vtkUnstructuredGrid>::New();
  d_grid_p->SetPoints(points);
}


void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<uint8_t> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (unsigned char i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<size_t> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (unsigned long i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<int> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (int i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<float> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (float i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(const std::string &name,
                                            const std::vector<double> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (double i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}

void rw::writer::VtkWriter::appendPointData(
    const std::string &name, const std::vector<util::Point3> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(3);
  array->SetName(name.c_str());

  array->SetComponentName(0, "x");
  array->SetComponentName(1, "y");
  array->SetComponentName(2, "z");

  double value[3];
  for (const auto &i : *data) {
    value[0] = i.d_x;
    value[1] = i.d_y;
    value[2] = i.d_z;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetPointData()->AddArray(array);
}



void rw::writer::VtkWriter::appendCellData(const std::string &name,
                                           const std::vector<float> *data) {
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetNumberOfComponents(1);
  array->SetName(name.c_str());

  double value[1];
  for (float i : *data) {
    value[0] = i;
    array->InsertNextTuple(value);
  }

  d_grid_p->GetCellData()->AddArray(array);
}



void rw::writer::VtkWriter::addTimeStep(const double &timestep) {
  auto t = vtkDoubleArray::New();
  t->SetName("TIME");
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, timestep);
  d_grid_p->GetFieldData()->AddArray(t);
}

void rw::writer::VtkWriter::close() {
  d_writer_p->SetInputData(d_grid_p);
  d_writer_p->SetDataModeToAppended();
  d_writer_p->EncodeAppendedDataOn();
  if (d_compressType == "zlib")
    d_writer_p->SetCompressorTypeToZLib();
  else
    d_writer_p->SetCompressor(0);
  d_writer_p->Write();
}

void rw::writer::VtkWriter::appendFieldData(const std::string &name,
                                            const double &data) {
  auto t = vtkDoubleArray::New();
  t->SetName(name.c_str());
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, data);
  d_grid_p->GetFieldData()->AddArray(t);
}

void rw::writer::VtkWriter::appendFieldData(const std::string &name,
                                            const float &data) {
  auto t = vtkDoubleArray::New();
  t->SetName(name.c_str());
  t->SetNumberOfTuples(1);
  t->SetTuple1(0, data);
  d_grid_p->GetFieldData()->AddArray(t);
}
