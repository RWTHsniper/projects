#ifndef POSTPROCESSOR_H_
#define POSTPROCESSOR_H_

#include "node.h"
#include "element.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkPropAssembly.h>
#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>

class postProcessor
{
    private:
        // VARIABLES
        double minT;       // min value of the Temperature field
        double maxT;       // max value of the Temperature field

        // METHODS
        void evaluateLimits(setting*, node*);
        void vtkVisualization(setting*, node*, element*, int iter);
    protected:

    public:
        void postProcessorControl(setting*, node*, element*, int iter);
};


#endif /* POSTPROCESSOR_H_ */
