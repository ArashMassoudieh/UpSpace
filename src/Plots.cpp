#include "Grid.h"

#ifdef Qt_version
#include <qapplication.h>
#include <qmainwindow.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qwt_plot.h>
#include <qdebug.h>
#include "heteval.h"
#include "Spectogram.h"
#endif // Qt_version

#include "vtk.h"
#include "Pathway.h"
#include "StringOP.h"

#ifdef QT_version
class PlotWindow : public QMainWindow
{
public:
	PlotWindow(CGrid *grid, QWidget *parent = NULL) :
		QMainWindow(parent)
	{
		plot = new Plot(grid);
		setCentralWidget(plot);

		QToolBar *toolBar = new QToolBar(this);

		QToolButton *btnLoad = new QToolButton(toolBar);

#if QT_VERSION >= 0x040000
		btnLoad->setText("Load SVG");
		btnLoad->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
		toolBar->addWidget(btnLoad);
#else
		btnLoad->setTextLabel("Load SVG");
		btnLoad->setUsesTextLabel(true);
#endif

		addToolBar(toolBar);

		//connect(btnLoad, SIGNAL(clicked()), plot, SLOT(loadSVG()));
	}
	Plot *plot;
};


void CGrid::show_K_field()
{
	PlotWindow *P = new PlotWindow(this);

	P->resize(600, 400);
	P->show();
}
#endif // QT_version

void CGrid::show_K_field_vtk(double z_factor)
{


}

/*void CGrid::showthings(vector<vtkSmartPointer<vtkActor>> actors, string filename)
{
	#ifdef QT_version
	main_window->get_ui()->ShowOutput->append("Creating renderer ...");
	#endif // QT_version
	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	for (int i=0; i<actors.size(); i++)
		renderer->AddActor(actors[i]);

	//renderWindow->SetFullScreen(true);
	renderWindow->SetSize(1600, 800);
	renderWindow->Render();

	if (filename != "")
	{
		#ifdef QT_version
		qDebug() << "Saving screen shot" << endl;
		main_window->get_ui()->ShowOutput->append("Saving Screen Shot ...");
		#endif // QT_version

		vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
			vtkSmartPointer<vtkWindowToImageFilter>::New();
		#ifdef QT_version
		qDebug() << "Setting Window Image Filter" << endl;
		#endif // QT_version
		windowToImageFilter->SetInput(renderWindow);
		windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
		windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
		windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
		windowToImageFilter->Update();

		vtkSmartPointer<vtkPNGWriter> writer =
			vtkSmartPointer<vtkPNGWriter>::New();
		writer->SetFileName(filename.c_str());
		writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		#ifdef QT_version
		qDebug() << "Writing the image"<<endl;
		#endif // QT_version
		writer->Write();
	}
	//->FullScreenOn();
	renderWindow->SetSize(1600, 800);
	renderWindow->Render();
	renderWindowInteractor->Start();

}*/

vtkSmartPointer<vtkActor> CGrid::get_K_field_vtk_pdt(double z_factor)
{
	vtkSmartPointer<vtkPoints> points_3 =
		vtkSmartPointer<vtkPoints>::New();

	double maxk = max_K();
	double xx, yy, zz;
	for (unsigned int x = 0; x < GP.nx; x++)
	{
		for (unsigned int y = 0; y < GP.ny; y++)
		{
			xx = x*GP.dx;
			yy = y*GP.dy;
			zz = p[x][y].K[0];
			points_3->InsertNextPoint(xx, yy, zz / maxk*z_factor);
		}
	}

	// Add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> inputPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	inputPolyData->SetPoints(points_3);

	// Triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay =
		vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
	delaunay->SetInput(inputPolyData);
#else
	delaunay->SetInputData(inputPolyData);
#endif
	delaunay->Update();
	vtkPolyData* outputPolyData = delaunay->GetOutput();

	double bounds[6];
	outputPolyData->GetBounds(bounds);

	// Find min and max z
	double minz = bounds[4];
	double maxz = bounds[5];

	std::cout << "minz: " << minz << std::endl;
	std::cout << "maxz: " << maxz << std::endl;

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		outputPolyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);
		//std::cout << "dcolor: "
		//	<< dcolor[0] << " "
		//	<< dcolor[1] << " "
		//	<< dcolor[2] << std::endl;
		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
		//std::cout << "color: "
		//	<< (int)color[0] << " "
		//	<< (int)color[1] << " "
		//	<< (int)color[2] << std::endl;

		colors_2->InsertNextTupleValue(color);
	}

	outputPolyData->GetPointData()->SetScalars(colors_2);


	//Append the two meshes
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	//appendFilter->AddInputData(polydata);
	//appendFilter->AddInputData(polydata_1);
	appendFilter->AddInputData(outputPolyData);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName("surface.vtp");
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();

	return actor;

}

void CGrid::write_K_field_to_vtp(string filename, double z_factor, bool _log)
{
	vtkSmartPointer<vtkPoints> points_3 =
		vtkSmartPointer<vtkPoints>::New();

	double maxk = max_K();
	double xx, yy, zz;
	vtkSmartPointer<vtkFloatArray> values =
		vtkSmartPointer<vtkFloatArray>::New();

	values->SetNumberOfComponents(1);
	if (_log)
		values->SetName("Hydraulic Conductivity");
	else
		values->SetName("Log Hydraulic Conductivity");

	for (unsigned int x = 0; x < GP.nx; x++)
	{
		for (unsigned int y = 0; y < GP.ny; y++)
		{
			xx = x*GP.dx;
			yy = y*GP.dy;
			zz = p[x][y].K[0];
			if (!_log)
			{
				float t[1] = { float(zz) };
				points_3->InsertNextPoint(xx, yy, zz / maxk*z_factor);
				values->InsertNextTupleValue(t);
			}
			else
			{
				float t[1] = { float(zz)};
				points_3->InsertNextPoint(xx, yy, log(*t/maxk)*z_factor);
				values->InsertNextTupleValue(t);
			}
		}
	}

	// Add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> inputPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	inputPolyData->SetPoints(points_3);

	// Triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay =
		vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
	delaunay->SetInput(inputPolyData);
#else
	delaunay->SetInputData(inputPolyData);
#endif
	delaunay->Update();
	vtkPolyData* outputPolyData = delaunay->GetOutput();

	double bounds[6];
	outputPolyData->GetBounds(bounds);

	// Find min and max z
	double minz = bounds[4];
	double maxz = bounds[5];

	std::cout << "minz: " << minz << std::endl;
	std::cout << "maxz: " << maxz << std::endl;

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		outputPolyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);
		//std::cout << "dcolor: "
		//	<< dcolor[0] << " "
		//	<< dcolor[1] << " "
		//	<< dcolor[2] << std::endl;
		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
		//std::cout << "color: "
		//	<< (int)color[0] << " "
		//	<< (int)color[1] << " "
		//	<< (int)color[2] << std::endl;

		colors_2->InsertNextTupleValue(color);
	}

	outputPolyData->GetPointData()->SetScalars(values);


	//Append the two meshes
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	//appendFilter->AddInputData(polydata);
	//appendFilter->AddInputData(polydata_1);
	appendFilter->AddInputData(outputPolyData);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();


}


vtkSmartPointer<vtkActor> CGrid::traj_vtk_pdt(int trajno, double z_factor, double offset)
{
	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	for (int i=0; i<Traj.paths[trajno].positions.size(); i++)
	{
		double p[3] = { Traj.paths[trajno].positions[i].x, Traj.paths[trajno].positions[i].y, Traj.paths[trajno].positions[i].v[0]*z_factor/max_v_x + offset };
		points->InsertNextPoint(p);
	}
	vtkSmartPointer<vtkPolyLine> polyLine =
		vtkSmartPointer<vtkPolyLine>::New();

	polyLine->GetPointIds()->SetNumberOfIds(Traj.paths[trajno].positions.size());
	for (unsigned int i = 0; i < Traj.paths[trajno].positions.size(); i++)
	{
		polyLine->GetPointIds()->SetId(i, i);
	}

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();

	// Add the points to the dataset
	polyData->SetPoints(points);

	// Add the lines to the dataset
	polyData->SetLines(cells);

	// Find min and max z
	double minz = min_v_x*z_factor/max_v_x + offset;
	double maxz = max_v_x*z_factor / max_v_x + offset;

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

	for (int i = 0; i < polyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		polyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);

		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}

		colors_2->InsertNextTupleValue(color);
	}

	polyData->GetPointData()->SetScalars(colors_2);

	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputData(polyData);
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	return actor;
}

vtkSmartPointer<vtkPolyData> CGrid::traj_vtk_pdt_vtp(int trajno, double z_factor, double offset, bool _log, bool _color)
{
	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkFloatArray> values =
		vtkSmartPointer<vtkFloatArray>::New();

	values->SetNumberOfComponents(1);
	values->SetName("vx");

	vtkSmartPointer<vtkFloatArray> vy =
		vtkSmartPointer<vtkFloatArray>::New();

	vy->SetNumberOfComponents(1);
	vy->SetName("vy");

	vtkSmartPointer<vtkFloatArray> tm =
		vtkSmartPointer<vtkFloatArray>::New();

	tm->SetNumberOfComponents(1);
	tm->SetName("arrival time");

	for (int i = 0; i<Traj.paths[trajno].positions.size(); i++)
	{

		if (!_log)
		{
			double t[1] = { float(Traj.paths[trajno].positions[i].v[0]) };
			double _vy[1] = { float(Traj.paths[trajno].positions[i].v[1]) };
			double _at[1] = { float(Traj.paths[trajno].positions[i].t) };
			double p[3] = { Traj.paths[trajno].positions[i].x, Traj.paths[trajno].positions[i].y, t[0] * z_factor / max_v_x + offset };
			points->InsertNextPoint(p);
			values->InsertNextTuple(t);
			vy->InsertNextTuple(_vy);
			tm->InsertNextTuple(_at);


		}
		else
		{
			double t[1] = { float(Traj.paths[trajno].positions[i].v[0]) };
			double _vy[1] = { float(Traj.paths[trajno].positions[i].v[1]) };
			double _at[1] = { float(Traj.paths[trajno].positions[i].t) };
			double p[3] = { Traj.paths[trajno].positions[i].x, Traj.paths[trajno].positions[i].y, float(log(fabs(Traj.paths[trajno].positions[i].v[0]) + 1e-12)) * z_factor / max_v_x + offset };
			points->InsertNextPoint(p);
			values->InsertNextTuple(t);
			vy->InsertNextTuple(_vy);
			tm->InsertNextTuple(_at);
		}
	}
	vtkSmartPointer<vtkPolyLine> polyLine =
		vtkSmartPointer<vtkPolyLine>::New();

	polyLine->GetPointIds()->SetNumberOfIds(Traj.paths[trajno].positions.size());
	for (unsigned int i = 0; i < Traj.paths[trajno].positions.size(); i++)
	{
		polyLine->GetPointIds()->SetId(i, i);
	}

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();

	// Add the points to the dataset
	polyData->SetPoints(points);

	// Add the lines to the dataset
	polyData->SetLines(cells);

	// Find min and max z
	double minz = min_v_x*z_factor / max_v_x + offset;
	double maxz = max_v_x*z_factor / max_v_x + offset;

	if (_log)
	{
		minz = log(fabs(min_v_x)+1e-12)*z_factor / max_v_x + offset;
		maxz = log(fabs(max_v_x) + 1e-12)*z_factor / max_v_x + offset;
	}

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	if (_color)
	{
		// Generate the colors for each point based on the color map
		vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
			vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors_2->SetNumberOfComponents(3);
		colors_2->SetName("Colors");

		for (int i = 0; i < polyData->GetNumberOfPoints(); i++)
		{
			double p[3];
			polyData->GetPoint(i, p);

			double dcolor[3];
			colorLookupTable->GetColor(p[2], dcolor);

			unsigned char color[3];
			for (unsigned int j = 0; j < 3; j++)
			{
				color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
			}

			colors_2->InsertNextTupleValue(color);
		}


	}

	polyData->GetPointData()->SetScalars(values);
	polyData->GetPointData()->AddArray(vy);
	polyData->GetPointData()->AddArray(tm);
	// Visualization

	return polyData;
}


vector<vtkSmartPointer<vtkActor>> CGrid::trajs_vtk_pdt(double z_factor, double offset)
{
	vector<vtkSmartPointer<vtkActor>> outactors;
	for (int i = 0; i < Traj.paths.size(); i++)
		outactors.push_back(traj_vtk_pdt(i, z_factor, offset));

	return outactors;
}

void CGrid::trajs_vtk_pdt_to_vtp(string filename, double z_factor, double offset, bool _log, bool _color)
{
	vector<vtkSmartPointer<vtkPolyData>> outputmappers;
	for (int i = 0; i < Traj.paths.size(); i++)
		outputmappers.push_back(traj_vtk_pdt_vtp(i, z_factor, offset, _log, _color));

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	for (int i=0; i<outputmappers.size(); i++)
		appendFilter->AddInputData(outputmappers[i]);
#endif
	appendFilter->Update();



	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	#ifdef QT_version
	main_window->get_ui()->ShowOutput->append("Writing vtp file... ");
	#endif // QT_version
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();

}





/*void CGrid::screen_shot(string filename)
{
	vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
		vtkSmartPointer<vtkWindowToImageFilter>::New();
	windowToImageFilter->SetInput(renderWindow);
	windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
	windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
	windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
	windowToImageFilter->Update();

	vtkSmartPointer<vtkPNGWriter> writer =
		vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	writer->Write();

	//renderWindow->Render();
	//renderer->ResetCamera();
	//renderWindow->Render();
	//renderWindowInteractor->Start();
}*/

/*void CGrid::screenshot_test()
{
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter(0.0, 0.0, 0.0);
		sphereSource->SetRadius(5.0);
		sphereSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(sphereSource->GetOutputPort());

		vtkSmartPointer<vtkActor> actor =
			vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);

		vtkSmartPointer<vtkRenderer> renderer =
			vtkSmartPointer<vtkRenderer>::New();
		vtkSmartPointer<vtkRenderWindow> renderWindow =
			vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->AddRenderer(renderer);
		renderWindow->SetAlphaBitPlanes(1); //enable usage of alpha channel

		vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
			vtkSmartPointer<vtkRenderWindowInteractor>::New();
		renderWindowInteractor->SetRenderWindow(renderWindow);

		renderer->AddActor(actor);
		renderer->SetBackground(1, 1, 1); // Background color white

		renderWindow->Render();

		// Screenshot
		vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
			vtkSmartPointer<vtkWindowToImageFilter>::New();
		windowToImageFilter->SetInput(renderWindow);
		windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
		windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
		windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
		windowToImageFilter->Update();

		vtkSmartPointer<vtkPNGWriter> writer =
			vtkSmartPointer<vtkPNGWriter>::New();
		writer->SetFileName("screenshot2.png");
		writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		writer->Write();

		renderWindow->Render();
		renderer->ResetCamera();
		renderWindow->Render();
		renderWindowInteractor->Start();


}*/

vtkSmartPointer<vtkPolyData> pathway_vtk_pdt_vtp(CPathway &pathway, double z_factor, double offset)
{
	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkFloatArray> values_vx =
		vtkSmartPointer<vtkFloatArray>::New();
	values_vx->SetNumberOfComponents(1);
	values_vx->SetName("vx");

	vtkSmartPointer<vtkFloatArray> values_u =
		vtkSmartPointer<vtkFloatArray>::New();
	values_u->SetNumberOfComponents(1);
	values_u->SetName("u");

	vtkSmartPointer<vtkFloatArray> values_z =
		vtkSmartPointer<vtkFloatArray>::New();
	values_z->SetNumberOfComponents(1);
	values_z->SetName("z");

	vtkSmartPointer<vtkFloatArray> values_t =
		vtkSmartPointer<vtkFloatArray>::New();
	values_t->SetNumberOfComponents(1);
	values_t->SetName("arrival_time");


	for (int i = 0; i<pathway.positions.size(); i++)
	{
			double vx[1] = { pathway.positions[i].v[0] };
			double u[1] = { pathway.positions[i].u };
			double z[1] = { pathway.positions[i].z };
			double p[3] = { pathway.positions[i].x , pathway.positions[i].y, pathway.positions[i].v[0]*z_factor + offset };
			points->InsertNextPoint(p);
			values_vx->InsertNextTuple(vx);
			values_z->InsertNextTuple(z);
			values_u->InsertNextTuple(u);
	}
	vtkSmartPointer<vtkPolyLine> polyLine =
		vtkSmartPointer<vtkPolyLine>::New();

	polyLine->GetPointIds()->SetNumberOfIds(pathway.positions.size());
	for (unsigned int i = 0; i < pathway.positions.size(); i++)
	{
		polyLine->GetPointIds()->SetId(i, i);
	}

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> cells =
		vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(polyLine);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();

	// Add the points to the dataset
	polyData->SetPoints(points);

	// Add the lines to the dataset
	polyData->SetLines(cells);

	polyData->GetPointData()->SetScalars(values_vx);
	polyData->GetPointData()->AddArray(values_u);
	polyData->GetPointData()->AddArray(values_z);
	polyData->GetPointData()->AddArray(values_t);

	// Visualization

	return polyData;
}

vtkSmartPointer<vtkPolyDataMapper> pathways_vtk_pdt_vtp(CPathwaySet *pathwayset, double z_factor, double offset)
{
	vector<vtkSmartPointer<vtkPolyData>> outarray;
	for (int i = 0; i < pathwayset->paths.size(); i++)
		outarray.push_back(pathway_vtk_pdt_vtp(pathwayset->paths[i], z_factor, offset));

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	for (int i = 0; i<outarray.size(); i++)
		appendFilter->AddInputData(outarray[i]);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	return mapper;

}

void write_vtk(vtkSmartPointer<vtkPolyDataMapper> mapper, string filename)
{
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();


}

void CGrid::write_K_solution_to_vtp(string filename, double z_factor, bool _log)
{
	vtkSmartPointer<vtkPoints> points_3 =
		vtkSmartPointer<vtkPoints>::New();

	double maxk = max_K();
	double xx, yy, zz;
	vtkSmartPointer<vtkFloatArray> values =
		vtkSmartPointer<vtkFloatArray>::New();

	vtkSmartPointer<vtkFloatArray> _H =
		vtkSmartPointer<vtkFloatArray>::New();

	vtkSmartPointer<vtkFloatArray> v =
		vtkSmartPointer<vtkFloatArray>::New();

	values->SetNumberOfComponents(1);
	_H->SetNumberOfComponents(1);
	v->SetNumberOfComponents(2);

	if (_log)
		values->SetName("Hydraulic Conductivity");
	else
		values->SetName("Log Hydraulic Conductivity");

	_H->SetName("Hydraulic Head");
	v->SetName("Velocity");

	for (unsigned int x = 0; x < GP.nx; x++)
	{
		for (unsigned int y = 0; y < GP.ny; y++)
		{
			xx = x*GP.dx;
			yy = y*GP.dy;
			zz = p[x][y].K[0];
			if (!_log)
			{
				float t[1] = { float(zz) };
				float vv[2] = { p[x][y].V[0], p[x][y].V[1] };
				float HH[1] = { H[x][y] };
				points_3->InsertNextPoint(xx, yy, zz / maxk*z_factor);
				values->InsertNextTupleValue(t);
				v->InsertNextTupleValue(vv);
				_H->InsertNextTupleValue(HH);

			}
			else
			{
				float t[1] = { float(zz) };
				float vv[2] = { p[x][y].V[0], p[x][y].V[1] };
				float HH[1] = { H[x][y] };
				points_3->InsertNextPoint(xx, yy, log(*t / maxk)*z_factor);
				values->InsertNextTupleValue(t);
				v->InsertNextTupleValue(vv);
				_H->InsertNextTupleValue(HH);

			}
		}
	}

	// Add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> inputPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	inputPolyData->SetPoints(points_3);

	// Triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay =
		vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
	delaunay->SetInput(inputPolyData);
#else
	delaunay->SetInputData(inputPolyData);
#endif
	delaunay->Update();
	vtkPolyData* outputPolyData = delaunay->GetOutput();

	double bounds[6];
	outputPolyData->GetBounds(bounds);

	// Find min and max z
	double minz = bounds[4];
	double maxz = bounds[5];

	std::cout << "minz: " << minz << std::endl;
	std::cout << "maxz: " << maxz << std::endl;

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		outputPolyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);
		//std::cout << "dcolor: "
		//	<< dcolor[0] << " "
		//	<< dcolor[1] << " "
		//	<< dcolor[2] << std::endl;
		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
//		std::cout << "color: "
//			<< (int)color[0] << " "
//			<< (int)color[1] << " "
//			<< (int)color[2] << std::endl;

		colors_2->InsertNextTupleValue(color);
	}

	outputPolyData->GetPointData()->SetScalars(values);
	outputPolyData->GetPointData()->AddArray(_H);
	outputPolyData->GetPointData()->AddArray(v);

	//Append the two meshes
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	//appendFilter->AddInputData(polydata);
	//appendFilter->AddInputData(polydata_1);
	appendFilter->AddInputData(outputPolyData);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();


}


void CGrid::write_C_to_vtp(string filename, double z_factor, bool _log, vector<double> t)
{
	for (int i = 0; i < t.size(); i++)
	{
		string filename_1 = split(filename, '.')[0] + "_" + numbertostring(i) +"." + split(filename, '.')[1];
		write_C_to_vtp(filename_1, z_factor, _log, t[i]);
	}
}

void CGrid::write_C_to_vtp(string filename, double z_factor, bool _log, double t)
{
	vtkSmartPointer<vtkPoints> points_3 =
		vtkSmartPointer<vtkPoints>::New();


	double xx, yy, zz;
	vtkSmartPointer<vtkFloatArray> values =
		vtkSmartPointer<vtkFloatArray>::New();


	values->SetNumberOfComponents(1);

	if (_log)
		values->SetName("Log Concentration");
	else
		values->SetName("Concentration");

	for (unsigned int x = 0; x < GP.nx; x++)
	{
		for (unsigned int y = 0; y < GP.ny; y++)
		{
			xx = x*GP.dx;
			yy = y*GP.dy;
			zz = p[x][y].C[int(t/dt)];
			if (!_log)
			{
				float CC[1] = { float(zz) };
				points_3->InsertNextPoint(xx, yy, zz * z_factor);
				values->InsertNextTupleValue(CC);
			}
			else
			{
				float CC[1] = { float(zz) };
				points_3->InsertNextPoint(xx, yy, log(max(zz,1e-25))*z_factor);
				values->InsertNextTupleValue(CC);
			}
		}
	}

	// Add the grid points to a polydata object
	vtkSmartPointer<vtkPolyData> inputPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	inputPolyData->SetPoints(points_3);

	// Triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay =
		vtkSmartPointer<vtkDelaunay2D>::New();
#if VTK_MAJOR_VERSION <= 5
	delaunay->SetInput(inputPolyData);
#else
	delaunay->SetInputData(inputPolyData);
#endif
	delaunay->Update();
	vtkPolyData* outputPolyData = delaunay->GetOutput();

	double bounds[6];
	outputPolyData->GetBounds(bounds);

	// Find min and max z
	double minz = bounds[4];
	double maxz = bounds[5];

	std::cout << "minz: " << minz << std::endl;
	std::cout << "maxz: " << maxz << std::endl;

	// Create the color map
	vtkSmartPointer<vtkLookupTable> colorLookupTable =
		vtkSmartPointer<vtkLookupTable>::New();
	colorLookupTable->SetTableRange(minz, maxz);
	colorLookupTable->Build();

	// Generate the colors for each point based on the color map
	vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors_2->SetNumberOfComponents(3);
	colors_2->SetName("Colors");

//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;

	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
	{
		double p[3];
		outputPolyData->GetPoint(i, p);

		double dcolor[3];
		colorLookupTable->GetColor(p[2], dcolor);
		//std::cout << "dcolor: "
		//	<< dcolor[0] << " "
		//	<< dcolor[1] << " "
		//	<< dcolor[2] << std::endl;
		unsigned char color[3];
		for (unsigned int j = 0; j < 3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
		//std::cout << "color: "
		//	<< (int)color[0] << " "
		//	<< (int)color[1] << " "
		//	<< (int)color[2] << std::endl;

		colors_2->InsertNextTupleValue(color);
	}

	outputPolyData->GetPointData()->SetScalars(values);

	//Append the two meshes
	vtkSmartPointer<vtkAppendPolyData> appendFilter =
		vtkSmartPointer<vtkAppendPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
	appendFilter->AddInputConnection(input1->GetProducerPort());
	appendFilter->AddInputConnection(input2->GetProducerPort());
#else
	//appendFilter->AddInputData(polydata);
	//appendFilter->AddInputData(polydata_1);
	appendFilter->AddInputData(outputPolyData);
#endif
	appendFilter->Update();


	// Visualization
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputConnection(appendFilter->GetOutputPort());
	//mapper->SetInputData(polydata_1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(5);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
		vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(mapper->GetInput());
	// This is set so we can see the data in a text editor.
	writer->SetDataModeToAscii();
	writer->Write();

}

void CGrid::clear()
{
	p.clear();
	vx = CMatrix();
	vy = CMatrix();
	H = CMatrix();
	C = CMatrix();
	Kv.matr.clear();
	KD.matr.clear();
	Kt.matr.clear();
	pset.paths.clear();
	vx_dist.clear();
	vy_dist.clear();
	sect_dist.clear();
}







