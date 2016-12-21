The MevisLab network performs
1. Cartilage segmentation
2. combines the output of bone segmentation to produce final output


Steps. 
1. Open the Original image in itkImageFileReader
2. Open Bone Segmentation output in itkImageFileReader2
2. Place seeds, adjust threshold (if required) and run RegionGrowingMacro for femoral cartilage segmentation
3. Place seeds, adjust threshold (if required) and run RegionGrowingMacro1 for tibial cartilage segmentation
4. View image in OrthoView2D
5. Specify output file name and location in itkImageFileWriter and Save file