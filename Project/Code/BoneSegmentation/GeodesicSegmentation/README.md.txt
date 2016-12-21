The executable is \build\Debug\GeodesicSegmentation.exe
Created using cmake and Visual Studio 2010 
The same can be done on user's machine


3 contiguous sections per bone type in the 3D volume (excluding ones without bone)
Section 1 - smaller area towards the beginning of the volume 
Section 2 - bigger bone area
Section 3 - smaller area towards the end of the volume

Slice Number counted from 1


To execute :

<Executable> <InputFile> <OutputFile> <Start Slice Number of Femur Section 1> <Seedx> <Seedy> <Start Slice Number of Tibia Section 1> <Seedx> <Seedy> <Start Slice Number of Femur Section 2> <End Slice Number of Femur Section 2> <Seedx> <Seedy> <Start Slice Number of Tibia Section 2> <End Slice Number of End Tibia 2> <Seedx> <Seedy> <End Slice Number of Femur Section 3> <Seedx> <Seedy> <End Slice Number of Tibia Section 3> <Seedx> <Seedy>


Example :

C:\Users\Gourab\Documents\MedicalImag\KneeMRISegmentation\BoneSegmentation\GeodesicSegmentation\build\Debug\GeodesicSegmentation.exe C:\User\Gourab\Documents\MedicalImag\data_train_fp\image-061.mhd C:\Users\Gourab\Documents\MedicalImag\Output\boneoutput-061.mhd 24 144 149 20 140 241 39 75 123 104 30 71 141 247 93 141 146 89 157 247




References Used :

Geodesic Active Contours Segmentation
http://www.itk.org/ItkSoftwareGuide.pdf (pg - 697)


Texture based Geodesic Active Contours (Lorigo)
http://imagine.enpc.fr/publications/papers/98miccai.pdf


ITK Sample code for Geodesic Active Contours Segmentation (used and modified for our MIA project)
http://www.itk.org/Doxygen/html/Examples_2Segmentation_2GeodesicActiveContourImageFilter_8cxx-example.html
http://www.apache.org/licenses/LICENSE-2.0.txt




