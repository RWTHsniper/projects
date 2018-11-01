#include "postProcessor.h"

/***************************************************************************************************
 * void postProcessor::postProcessorControl(setting*, node*, element*)
 ***************************************************************************************************
 * This routine controls the post processing operations.
 **************************************************************************************************/
void postProcessor::postProcessorControl(setting* stngs, node* nodes, element* elems, int iter)
{
    cout << endl;
    cout << "====== Post-processing ======================================================" << endl;

    evaluateLimits(stngs, nodes);
    vtkVisualization(stngs, nodes, elems, iter);

    return;
}


/***************************************************************************************************
 * void postProcessor::evaluateLimits(setting*, node*)
 ***************************************************************************************************
 * Evaluates the maximum and minimum temperatures in the field
 **************************************************************************************************/
void postProcessor::evaluateLimits(setting* stngs, node* nodes)
{
    int nn = stngs->getNn();    // Number of nodes
    double T;   // Temporary temperature value
    minT = std::numeric_limits<double>::max();
    maxT = std::numeric_limits<double>::min();

    for (int i = 0; i < nn; i++)
    {
        T = nodes->getT(i);
            if(T < minT)
                minT = T;
            if(T > maxT)
                maxT = T;
    }
    cout << "Tmin = " << minT << endl;
    cout << "Tmax = " << maxT << endl;

    return;
}


template <typename T>
std::string to_string(T value)
{
	std::ostringstream os ;
	os << value ;
	return os.str() ;
}


/***************************************************************************************************
 * void postProcessor::vtkVisualization(setting*, node*, element*)
 ***************************************************************************************************
 * Main visualization function
 **************************************************************************************************/
void postProcessor::vtkVisualization(setting* stngs, node* nodes, element* elems, int iter)
{// Visualizing connectivity, temperature
    int nn  = stngs->getNn();
    int ne  = stngs->getNe();
    int nen = elems->getNen();
    string dummy;

    // VTK Double Array
    vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
    pcoords->SetNumberOfComponents(nen);
    pcoords->SetNumberOfTuples(nn);

    ///vtkDoubleArray type pcoords is filled with the data in meshPoints.
    for (int i=0; i<nn; i++)
        pcoords->SetTuple3(i,nodes->getX(i),nodes->getY(i),0.0f);

    ///vtkPoints type outputPoints is filled with the data in pcoords.
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetData(pcoords);

    ///Connectivity is written to vtkCellArray type outputCells
    vtkSmartPointer<vtkCellArray> connectivity = vtkSmartPointer<vtkCellArray>::New();
    for(int i=0; i<ne; i++)
    {
        connectivity->InsertNextCell(nen);
        for(int j=0; j<nen; j++)
                connectivity->InsertCellPoint(elems->getConn(i,j));
    }

    // Scalar property
    vtkSmartPointer<vtkDoubleArray> temperature = vtkSmartPointer<vtkDoubleArray>::New();
    temperature->SetName("Temperature");
    for(int i=0; i<nn; i++)
        temperature->InsertNextValue(nodes->getT(i));


    // Previously collected data which are outputPoints, outputCells, scalarProperty, are written to
    // vtkPolyData type polydata var.
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(outputPoints);
    polydata->SetPolys(connectivity);
    polydata->GetPointData()->SetScalars(temperature);

    // Whatever collected in polydata above is written to the vtk file below.
    // vtkDataSetWriter is for leagacy VTK format, vtkXMLDataSetWriter is for VTK XML format.
    // vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    dummy = stngs->getTitle();
    dummy.append(".");
    dummy.append(to_string(iter));  iter++; cout << "inside postP iter: " << iter << endl;
    dummy.append(".vtk");
    writer->SetFileName(dummy.c_str());
    writer->SetInputData(polydata);
    writer->Write();

    // ///In the below section, data collected in polydata var is rendered.
    // // Create the mapper and set the appropriate scalar range. (default is (0,1)
    // vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    // mapper->SetInput(polydata);
    // mapper->SetScalarRange(minT, maxT);

    // // Create an actor.
    // vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    // actor->SetMapper(mapper);

    // // Scalar bar
    // vtkSmartPointer<vtkScalarBarActor> scalarBar =
    // vtkSmartPointer<vtkScalarBarActor>::New();
    // scalarBar->SetLookupTable(mapper->GetLookupTable());
    // scalarBar->SetTitle("Temperature");
    // scalarBar->SetNumberOfLabels(10);

    // // Create a lookup table to share between the mapper and the scalarbar
    // vtkSmartPointer<vtkLookupTable> hueLut =
    // vtkSmartPointer<vtkLookupTable>::New();
    // hueLut->SetTableRange (0, 1);
    // hueLut->SetHueRange (0.7, 0.0);
    // hueLut->SetSaturationRange (1, 1);
    // hueLut->SetValueRange (1, 1);
    // hueLut->Build();

    // // llokup table shared between mapper and scalarbar
    // mapper->SetLookupTable( hueLut );
    // scalarBar->SetLookupTable( hueLut );

    // // Create the rendering objects.
    // vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    // renderer->GradientBackgroundOn();
    // renderer->SetBackground(1,1,1);
    // renderer->SetBackground2(0,0,0);

    // vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    // renWin->AddRenderer(renderer);
    // renWin->SetSize(1000,600);

    // cout << "> Rendering window is initiated." << endl;
    // vtkSmartPointer<vtkRenderWindowInteractor> iren = 
    // vtkSmartPointer<vtkRenderWindowInteractor>::New();
    // iren->SetRenderWindow(renWin);

    // vtkSmartPointer<vtkAxesActor> axes =
    // vtkSmartPointer<vtkAxesActor>::New();

    // vtkSmartPointer<vtkOrientationMarkerWidget> widget =
    // vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    // widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    // widget->SetOrientationMarker( axes );
    // widget->SetInteractor( iren );
    // widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    // widget->SetEnabled( 1 );
    // widget->InteractiveOff();

    // renderer->AddActor(actor);
    // renderer->AddActor2D(scalarBar);

    // renderer->ResetCamera();
    // renWin->Render();

    // // Begin mouse interaction
    // iren->Start();

    return;
}

