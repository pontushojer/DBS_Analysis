#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <iostream>
#include <string>

using namespace cv;
using namespace std;

static void help()
{
    cout << "\nThis program demonstrates circle finding with the Hough transform.\n"
            "Usage:\n"
            "./houghcircles <input_image_name> <output_image_name>\n" << endl
            << "If not specified the output will be out.<input_image_name>" << endl;
}

int main(int argc, char** argv)
{
    if(argc < 2 or argc > 3) 
    {
        help();
        return -1;
    }
    
    const string in_filename = argv[1];
    string out_filename = "out."+in_filename;
    
    if (argc == 3)
        out_filename = argv[2];
        

    Mat img = imread(in_filename, 0);
    if(img.empty())
    {
        help();
        cout << "can not open " << in_filename << endl;
        return -1;
    }

    Mat cimg;
    medianBlur(img, img, 5);        // amount of blurring to compensate for non perfect circles. PLAY WITH THIS!
    cvtColor(img, cimg, COLOR_GRAY2BGR);

    vector<Vec3f> circles;
    HoughCircles(img, circles, HOUGH_GRADIENT,
                 1,     // inverse ratio of resolution
                 10,    // minimum distance between centers     PLAY WITH THIS!
                 200,   // upper threshold for the internal Canny edge detector
                 25,    // threshold for center detection       
                 1,     // minimum radius           PLAY WITH THIS!
                 30     // maximum radius           PLAY WITH THIS!
                 );
    
    double total = 0;
   
    for( size_t i = 0; i < circles.size(); i++ )
    {
        Vec3i c = circles[i];
        circle( cimg, Point(c[0], c[1]), c[2], Scalar(0,0,255), 1, LINE_AA);        // circumference: ... scalar(blue,green,red), line weight
        circle( cimg, Point(c[0], c[1]), 2, Scalar(0,255,0), 3, LINE_AA);           // center point
        std::cout << "Found circle at coordinates " << c[0] << "." << c[1] << " of radius " << c[2] << " pixels." << std::endl;
        total += c[2];
    }
    
    double average = total/circles.size();
        
    std::cout << "\nImage sample: "+in_filename << std::endl;
    std::cout << "Droplets found: " << circles.size() << std::endl;
    std::cout << "Average droplet size in pixels: " << average << std::endl;
    std::cout << "Average droplet size in um: " << average*1.12352 << std::endl;
    
    imshow("detected circles", cimg);
    imwrite(out_filename, cimg);
    // while(waitKey() != 'q') {}       // close image and terminate only when pressing q (use instead of next line).
    //waitKey();                          // wait to close image and terminate script until any key is pressed 

    return 0;
}
