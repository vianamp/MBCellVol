// ==============================================================
// MBCellVol: This routine calculates the cellular volume based
// the mitochondrial content. The hypothesys is that the boundary
// of mitochondrial network is a good approximation for the cell
// boundary.
// ==============================================================

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vtkMath.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkLongArray.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkTransform.h>
#include <vtkDelaunay3D.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkVectorText.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkKdTree.h>
#include <vtkPoints.h>

bool _export_nodes_projection = false;
bool _export_ellipsoid_polydata = true;

void SavePolyData(vtkPolyData *PolyData, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving PolyData from XYZ list...\n");
    #endif

    #ifdef DEBUG
        printf("\t#Points in PolyData file: %llu.\n",(vtkIdType)PolyData->GetNumberOfPoints());
    #endif

    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer -> SetFileType(VTK_BINARY);
    Writer -> SetFileName(FileName);
    Writer -> SetInputData(PolyData);
    Writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif
}

template<typename sstype> void swap(sstype *x, sstype *y) {
     sstype t = *y; *y = *x; *x = t;
}

void sort(double v1, double v2, double v3, int *i1, int *i2, int *i3) {
    double l1 = v1, l2 = v2, l3 = v3;
    if (fabs(l1) > fabs(l2)) { swap(&l1,&l2); swap(i1,i2); }
    if (fabs(l2) > fabs(l3)) { swap(&l2,&l3); swap(i2,i3); }
    if (fabs(l1) > fabs(l2)) { swap(&l1,&l2); swap(i1,i2); }
}

void CalculateNodes2DProjection(vtkPolyData *Ellipsoid, const char FileName[]) {
    #ifdef DEBUG
        printf("Calculating 2D projection...\n");
    #endif

    int j, N = 50;
    double d, r[3];
    float x, y, z, fx, fy;
    char _fullpath[256];
    sprintf(_fullpath,"%s.coo",FileName);

    vtkPoints *NP = vtkPoints::New();

    vtkPoints *Points = Ellipsoid -> GetPoints();
    vtkSmartPointer<vtkKdTree> KdTree = vtkSmartPointer<vtkKdTree>::New();
    KdTree -> BuildLocatorFromPoints(Points);

    FILE *fcoo = fopen(_fullpath,"r");

    sprintf(_fullpath,"%s_2dproj.coo",FileName);
    FILE *fcoo2d = fopen(_fullpath,"w");
    while (fscanf(fcoo,"%f %f %f",&x,&y,&z) != EOF) {
        r[0] = 0.056*x; r[1] = 0.056*y; r[2] = 0.2*z;
        j = KdTree -> FindClosestPoint(r,d);
        Points -> GetPoint(j,r);

        x = ((j-2)%((N-2)*N))/(N-2);
        y = (j-2)%(N-2);
        // Add small random factor to avoid nodes with same coordinates
        fx = double(x)/(N-1);
        fy = double(y)/(N-3);

        r[0] = x; r[1] = y; r[2] = 0.0;

        NP -> InsertNextPoint(r);

        fprintf(fcoo2d,"%1.3f\t%1.3f\t0.0\n",fx,fy);
    }
    fclose(fcoo);
    fclose(fcoo2d);

    vtkPolyData *NPoly = vtkPolyData::New();
    NPoly -> SetPoints(NP);

    SavePolyData(NPoly,"ProjectedNodes.vtk");
 
    #ifdef DEBUG
        printf("\tDone!\n");
    #endif
}

void GetEllipsoidFrom3DConvexHull(const char FileName[], double Rad[3]) {

    char _fullpath[256];
    sprintf(_fullpath,"%s_surface.vtk",FileName);

    #ifdef DEBUG
        printf("Loading surface from %s\n",_fullpath);
    #endif

    vtkSmartPointer<vtkPolyDataReader> PolyReader = vtkSmartPointer<vtkPolyDataReader>::New();    
    PolyReader -> SetFileName(_fullpath);
    PolyReader -> Update();

    vtkPolyData *Surface = PolyReader -> GetOutput();

    #ifdef DEBUG
        printf("Calculating Delaunay triangulation...\n");
    #endif

    vtkSmartPointer<vtkDelaunay3D> Delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();
    Delaunay3D -> GlobalWarningDisplayOff();
    Delaunay3D -> SetAlpha(0.0);
    Delaunay3D -> SetTolerance(0.001);
    Delaunay3D -> SetOffset(2.5);
    Delaunay3D -> SetInputData(Surface);
    Delaunay3D -> Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> UnStrctToPolyData = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    UnStrctToPolyData -> SetInputData(Delaunay3D->GetOutput());
    UnStrctToPolyData -> Update();

    vtkPolyData *SmoothCVXH = UnStrctToPolyData -> GetOutput();

    vtkPoints *CVXHPoints = SmoothCVXH -> GetPoints();

    #ifdef DEBUG
        printf("Fitting ellipsoid using %lld points...\n",CVXHPoints->GetNumberOfPoints());
    #endif

    /*Fitting  the set of  points  CVXHPoints by  an arbitrary ellipsoid.
    This code is  based on  the script by Yuri  Petrov for  Matlab, which
    can be found  at
    http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
    Given  a  set  of  points [X,Y,Z], we try to fit  them by the surface
    described by equation

        ax^2 + by^2 + cz^2 + 2dxy + 2exz + 2fyz + 2gx + 2hy + 2iz = 1

    The coefficients are found by solving the linear system
                                A'Ax = A'B,
    where
                                A = [x^2 y^2 z^2 2xy 2xz 2yz 2x 2y 2z].

    On this code, ATA/ATB denotes A'A/A'B.

    --Matheus Palhares Viana - 8.30.2013 - vianamp@gmail.com*/

    int i, j, k;
    int n = CVXHPoints -> GetNumberOfPoints();

    double* ATB = new double[9];
    double**ATA = new double*[9];
    for (i = 0; i < 9; i++) ATA[i] = new double[9];

    double **A = new double*[n];
    for (i = 0; i < n; i++) A[i] = new double[9];

    for (i = 0; i < 9; i++) {
        ATB[i] = 0.0;
        for (j = 0; j < 9; j++) ATA[i][j] = 0.0;
    }

    double x, y, z, r[3];
    for (i = 0; i < n; i++) {
        CVXHPoints -> GetPoint(i,r);
        x = r[0]; y=r[1]; z = r[2];
        A[i][0] =     x*x; A[i][1] =     y*y; A[i][2] =     z*z;
        A[i][3] = 2.0*x*y; A[i][4] = 2.0*x*z; A[i][5] = 2.0*y*z;
        A[i][6] =   2.0*x; A[i][7] =   2.0*y; A[i][8] =   2.0*z;
    }

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            for (k = 0; k < n; k++) ATA[i][j] += A[k][i] * A[k][j];
        }
    }

    for (i = 0; i < 9; i++) {
        for (k = 0; k < n; k++) ATB[i] += A[k][i];
    }

    #ifdef DEBUG
        printf("\tSolving the 1st linear system...\n");
    #endif

    vtkMath::SolveLinearSystem(ATA,ATB,9);

    double **C = new double*[4];
    for (i = 0; i < 4; i++) C[i] = new double[4];

    C[0][0] = ATB[0]; C[0][1] = ATB[3]; C[0][2] = ATB[4]; C[0][3] = ATB[6];
    C[1][0] = ATB[3]; C[1][1] = ATB[1]; C[1][2] = ATB[5]; C[1][3] = ATB[7];
    C[2][0] = ATB[4]; C[2][1] = ATB[5]; C[2][2] = ATB[2]; C[2][3] = ATB[8];
    C[3][0] = ATB[6]; C[3][1] = ATB[7]; C[3][2] = ATB[8]; C[3][3] = -1;

    double* B2 = new double[3];
    double**A2 = new double*[3];
    for (i = 0; i < 3; i++) A2[i] = new double[3];

    for (i = 0; i < 3; i++) {
        B2[i] = ATB[6+i];
        for (j = 0; j < 3; j++) A2[i][j] = -C[i][j];
    }

    #ifdef DEBUG
        printf("\tSolving the 2nd linear system...\n");
    #endif

    vtkMath::SolveLinearSystem(A2,B2,3);

    double EllipsoidCenter[] = {B2[0],B2[1],B2[2]};

    double **T = new double*[4];
    double **R = new double*[4];
    double **Q = new double*[4];
    for (i = 0; i < 4; i++) {
        T[i] = new double[4];
        R[i] = new double[4];
        Q[i] = new double[4];
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            T[i][j] = R[i][j] = Q[i][j] = 0.0;
        }
    }

    for (i = 0; i < 4; i++) T[i][i] = 1.0;

    T[3][0] = EllipsoidCenter[0]; T[3][1] = EllipsoidCenter[1]; T[3][2] = EllipsoidCenter[2];

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) Q[i][j] += C[i][k] * T[j][k];
        }
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) R[i][j] += T[i][k] * Q[k][j];
        }
    }

    double V[3][3], EVE[3][3], EVA[3];

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            V[i][j] = -R[i][j]/R[3][3];
        }
    }

    #ifdef DEBUG
        printf("\tDiagonalizing the eigensystem...\n");
    #endif

    vtkMath::Diagonalize3x3(V,EVA,EVE);

    int l1 = 0, l2 = 1, l3 = 2;
    sort(EVA[0],EVA[1],EVA[2],&l1,&l2,&l3);

    #ifdef DEBUG
        printf("\t\tEigenvalues: %1.3f\t%1.3f\t%1.3f\n",EVA[l1],EVA[l2],EVA[l3]);    
    #endif

    // Ellipsoid axes
    int npixelstoadd = 0;
    Rad[0] = npixelstoadd+sqrt(1.0/EVA[l3]);
    Rad[1] = npixelstoadd+sqrt(1.0/EVA[l2]);
    Rad[2] = npixelstoadd+sqrt(1.0/EVA[l1]);

    if (_export_ellipsoid_polydata) {

        double vx[] = {1,0,0};
        double vy[] = {0,1,0};
        double vz[] = {0,0,1};
        double ex[] = {EVE[0][l3],EVE[1][l3],EVE[2][l3]};
        double ey[] = {EVE[0][l2],EVE[1][l2],EVE[2][l2]};
        double ez[] = {EVE[0][l1],EVE[1][l1],EVE[2][l1]};
        double q_vx_ex = acos(ex[0]);

        vtkSmartPointer<vtkSphereSource> Sphere = vtkSmartPointer<vtkSphereSource>::New();
        Sphere -> SetRadius(Rad[2]);
        Sphere -> SetPhiResolution(50);
        Sphere -> SetThetaResolution(50);
        Sphere -> Update();

        vtkSmartPointer<vtkTransform> Scaling = vtkSmartPointer<vtkTransform>::New();
        Scaling -> Scale(Rad[0]/Rad[2],Rad[1]/Rad[2],1.0);
        Scaling -> Update();

        vtkSmartPointer<vtkTransformPolyDataFilter> ImplementScalingT = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        ImplementScalingT -> SetInputData(Sphere->GetOutput());
        ImplementScalingT -> SetTransform(Scaling);
        ImplementScalingT -> Update();

        #ifdef DEBUG
            printf("\tGenerating ellipsoid...\n");
        #endif

        vtkPolyData *Ellipsoid = ImplementScalingT -> GetOutput();

        #ifdef DEBUG
            printf("\tTranslating the ellipsoid to the correct position...\n");
        #endif

        /*Moving the ellipsoid to the right place in 3D space.This is based on
        two  rotations. The first  one moves  the x versor  towards the  vector
        ex, which is  the ellipsoid major axis.  This  rotation must be done
        around the  ortogonal versor  between  x and ex. The  second rotation
        moves  the now  rotated  secondmajor  axis towards its right position.
        This rotation must be done around the ellipsoid major axis.
        */

        double ox[3];
        vtkMath::Cross(vx,ex,ox);
        vtkMath::Normalize(ox);
        double nn = sqrt(ox[1]*ox[1]+ox[2]*ox[2]);
        vtkSmartPointer<vtkTransform> RotOp = vtkSmartPointer<vtkTransform>::New();
        RotOp -> RotateWXYZ(180.0/3.1415*q_vx_ex,ox);
        RotOp -> Update();

        vtkSmartPointer<vtkTransformPolyDataFilter> ImplementRotOpT = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        ImplementRotOpT -> SetInputData(Ellipsoid);
        ImplementRotOpT -> SetTransform(RotOp);
        ImplementRotOpT -> Update();
        Ellipsoid = ImplementRotOpT -> GetOutput();

        /* Position of ey after rotation RotOp. This equation can be found
        at http://mathworld.wolfram.com/RodriguesRotationFormula.html.
        (Rodrigues rotation matrix) */

        x=ox[0];y=ox[1];z=ox[2];
        double RD[3][3] = {{cos(q_vx_ex)+x*x*(1-cos(q_vx_ex)),x*y*(1-cos(q_vx_ex))-z*sin(q_vx_ex),x*z*(1-cos(q_vx_ex))+y*sin(q_vx_ex)},
                           {y*x*(1-cos(q_vx_ex))+z*sin(q_vx_ex),cos(q_vx_ex)+y*y*(1-cos(q_vx_ex)),y*z*(1-cos(q_vx_ex))-x*sin(q_vx_ex)},
                           {z*x*(1-cos(q_vx_ex))-y*sin(q_vx_ex),z*y*(1-cos(q_vx_ex))+x*sin(q_vx_ex),cos(q_vx_ex)+z*z*(1-cos(q_vx_ex))}};

        double rotvy[3] = {RD[0][1]*vy[1],RD[1][1]*vy[1],RD[2][1]*vy[1]};

        double q_ey_rotvy = -acos(vtkMath::Dot(rotvy,ey));

        /* Checking if the second rotation must be done clockwise or not */
        x=ex[0];y=ex[1];z=ex[2];
        double RD2[3][3] = {{cos(q_ey_rotvy)+x*x*(1-cos(q_ey_rotvy)),x*y*(1-cos(q_ey_rotvy))-z*sin(q_ey_rotvy),x*z*(1-cos(q_ey_rotvy))+y*sin(q_ey_rotvy)},
                    {y*x*(1-cos(q_ey_rotvy))+z*sin(q_ey_rotvy),cos(q_ey_rotvy)+y*y*(1-cos(q_ey_rotvy)),y*z*(1-cos(q_ey_rotvy))-x*sin(q_ey_rotvy)},
                    {z*x*(1-cos(q_ey_rotvy))-y*sin(q_ey_rotvy),z*y*(1-cos(q_ey_rotvy))+x*sin(q_ey_rotvy),cos(q_ey_rotvy)+z*z*(1-cos(q_ey_rotvy))}};

        double rotrotvy[3] = {0,0,0};
        for (i=0;i<3;i++) {
            rotrotvy[0] += RD2[0][i]*rotvy[i];
            rotrotvy[1] += RD2[1][i]*rotvy[i];
            rotrotvy[2] += RD2[2][i]*rotvy[i];
        }

        double q_rrvy_ey = vtkMath::Dot(rotrotvy,ey);
        if ((1.0-q_rrvy_ey)>1E-5) {
            q_ey_rotvy *= -1.0;
        }

        /* Performing the second Rotation */
        vtkSmartPointer<vtkTransform> RotOp2 = vtkSmartPointer<vtkTransform>::New();
        RotOp2 -> RotateWXYZ(180.0/3.1415*q_ey_rotvy,ex);
        RotOp2 -> Update();

        vtkSmartPointer<vtkTransformPolyDataFilter> ImplementRotOp2T = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        ImplementRotOp2T -> SetInputData(Ellipsoid);
        ImplementRotOp2T -> SetTransform(RotOp2);
        ImplementRotOp2T -> Update();
        Ellipsoid = ImplementRotOp2T -> GetOutput();

        vtkSmartPointer<vtkTransform> TransOp = vtkSmartPointer<vtkTransform>::New();
        TransOp -> Translate(EllipsoidCenter[0],EllipsoidCenter[1],EllipsoidCenter[2]);
        TransOp -> Update();

        vtkSmartPointer<vtkTransformPolyDataFilter> ImplementTransOp = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        ImplementTransOp -> SetInputData(Ellipsoid);
        ImplementTransOp -> SetTransform(TransOp);
        ImplementTransOp -> Update();
        Ellipsoid = ImplementTransOp -> GetOutput();

        sprintf(_fullpath,"%s_ellipsoid.vtk",FileName);
        SavePolyData(Ellipsoid,_fullpath);

        vtkSmartPointer<vtkPoints> EVE_Points = vtkSmartPointer<vtkPoints>::New();
        EVE_Points -> InsertPoint(0,EllipsoidCenter[0],EllipsoidCenter[1],EllipsoidCenter[2]);
        EVE_Points -> InsertPoint(1,EllipsoidCenter[0]+Rad[0]*ex[0],EllipsoidCenter[1]+Rad[0]*ex[1],EllipsoidCenter[2]+Rad[0]*ex[2]);
        EVE_Points -> InsertPoint(2,EllipsoidCenter[0]+Rad[1]*ey[0],EllipsoidCenter[1]+Rad[1]*ey[1],EllipsoidCenter[2]+Rad[1]*ey[2]);
        EVE_Points -> InsertPoint(3,EllipsoidCenter[0]+Rad[2]*ez[0],EllipsoidCenter[1]+Rad[2]*ez[1],EllipsoidCenter[2]+Rad[2]*ez[2]);
        EVE_Points -> InsertPoint(4,EllipsoidCenter[0]+Rad[2]*rotvy[0],EllipsoidCenter[1]+Rad[2]*rotvy[1],EllipsoidCenter[2]+Rad[2]*rotvy[2]);
        vtkSmartPointer<vtkCellArray> EVE_CellArray = vtkSmartPointer<vtkCellArray>::New();
        EVE_CellArray -> InsertNextCell(2);
        EVE_CellArray -> InsertCellPoint(0); EVE_CellArray -> InsertCellPoint(1);
        EVE_CellArray -> InsertNextCell(2);
        EVE_CellArray -> InsertCellPoint(0); EVE_CellArray -> InsertCellPoint(2);
        EVE_CellArray -> InsertNextCell(2);
        EVE_CellArray -> InsertCellPoint(0); EVE_CellArray -> InsertCellPoint(3);
        EVE_CellArray -> InsertNextCell(2);
        EVE_CellArray -> InsertCellPoint(0); EVE_CellArray -> InsertCellPoint(4);
        vtkSmartPointer<vtkPolyData> EVEVectors = vtkSmartPointer<vtkPolyData>::New();
        EVEVectors -> SetPoints(EVE_Points);
        EVEVectors -> SetLines(EVE_CellArray);

        sprintf(_fullpath,"%s_ellipsoid_axes.vtk",FileName);
        SavePolyData(EVEVectors,_fullpath);

        if (_export_nodes_projection) {
            CalculateNodes2DProjection(Ellipsoid,FileName);
        }

    }

    /* Clean memory */

    for (i = 0; i < n; i++) {
        delete[] A[i];
        if (i<9) {
            delete[] ATA[i];
            if (i<4) {
                delete[] C[i];
                delete[] T[i];
                delete[] R[i];
                delete[] Q[i];
                if (i<3) {
                    delete[] A2[i];
                }
            }
        }
    }
    delete[] A;
    delete[] ATB;
    delete[] ATA;
    delete[] C;
    delete[] T;
    delete[] R;
    delete[] Q;
    delete[] B2;
    delete[] A2;

}


/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     

    int i;
    char _impath[128];
    sprintf(_impath,"");

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-export_off")) {
            _export_ellipsoid_polydata = false;
        }
        if (!strcmp(argv[i],"-project_nodes")) {
            _export_nodes_projection = true;
        }
    }

    // Generating summary file and writing the header
    char _summaryfilename[256];
    sprintf(_summaryfilename,"%scell_volume.txt",_impath);
    FILE *fsummary = fopen(_summaryfilename,"w");
    fprintf(fsummary,"MBCellVol V1.0\n");
    fprintf(fsummary,"Folder: %s\n",_impath);
    fprintf(fsummary,"Cell path\tMinor radius (um)\tMid radius (um)\tMajor axis (um)\tSurface Area (um2)\tVolume (um3)\n");
    fclose(fsummary);

    double Rad[3], SA, V;
    char _tifffilename[256];
    char _tifflistpath[128];
    sprintf(_tifflistpath,"%smitograph.files",_impath);
    FILE *f = fopen(_tifflistpath,"r");
    while (fgets(_tifffilename,256, f) != NULL) {

        _tifffilename[strcspn(_tifffilename, "\n" )] = '\0';
        GetEllipsoidFrom3DConvexHull(_tifffilename,Rad);

        SA = 4*3.141592*pow( pow(Rad[2]*Rad[1],1.6) + pow(Rad[2]*Rad[0],1.6) + pow(Rad[1]*Rad[0],1.6) ,1.0/1.6);

        V = 4.0/3.0 * 3.141592 * Rad[0]*Rad[1]*Rad[2];

        fsummary = fopen(_summaryfilename,"a");
        fprintf(fsummary,"%s\t%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\n",_tifffilename,Rad[0],Rad[1],Rad[2],SA,V);
        fclose(fsummary);

        printf("%s\n",_tifffilename);
      
    }

}
