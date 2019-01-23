#ifdef QT_version
#include "heteval.h"
#include <QtWidgets/QApplication>
#include <sys/resource.h>
#endif // Qt_version
#include <iostream>
#include "Grid.h"
#include "2DMap.h"
#include "NormalDist.h"

int main(int argc, char *argv[])
{

	/*CTimeSeries X;
	for (int i=0; i<100000; i++)
        X.append(getnormalrand(5,3));

    CTimeSeries X_dist = X.distribution_fw(400);
    CTimeSeries X_dist_smooth = X.distribution_fw(400,0.2);
    X_dist.writefile("X_dist.txt");
    X_dist_smooth.writefile("X_dist_s.txt");*/

    /*TDMap M(10,0,1);
    for (int k=0; k<100000; k++)
    {
        unsigned int i=unitrandom()*10;
        unsigned int j=unitrandom()*10;
        M.add_val(i,j,1);
    }

    M.normalize();
    M.writetofile_GNU("map.gnu");*/

    /*CTimeSeries X;
    for (int i=0; i<100000; i++)
    {
        X.append(exp(getnormalrand(0,1)));
    }

    CTimeSeries dist = X.distribution(40, 0);
    dist.writefile("dist.txt");
	CTimeSeries unif_cum = X.getcummulative_direct(100);
	unif_cum.writefile("unif_cum.txt");
	X.rank_bd(100).writefile("rank_bd.txt");*/

	#ifdef QT_version
    QApplication a(argc, argv);
    HETEVAL w;
    w.showMaximized();
    return a.exec();
    #else
    string filename;
    cout<<"Enter the input file name: ";
    cin>>filename;
    //cout << "reading [" << filename << "]..." << endl;
    //filename = "/home/arash/Projects/Upscaling_outputs/input_test_BTC_log_normal_loop_std1_l_res_corr_s.txt";
    //filename = "/home/arash/Projects/UpscalingInputfiles/input_test_BTC_log_normal_loop_std1_l_res_corr.txt";
///xxx///
    filename = "/home/arash/Projects/UpscalingInputs/test.txt";
    CGrid G(filename);
    cout << "running [" << filename << "]..." << endl;
    G.runcommands_qt();
    return 0;
    #endif // Qt_version
}

//filename = "/home/arash/Projects/Upscaling_outputs/input_test_BTC_log_normal_loop_std1_l_iso_res_corr_s_1.txt";
