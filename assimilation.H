		// find ID of patch
		label patchID = mesh.boundaryMesh().findPatchID(interfaceName);
		// check patch has been found
		if(patchID == -1)
		{
				FatalError << "patch not found!" << exit(FatalError);
		}

		//calculate the heat flux of the inerface
		volScalarField& he = thermo.he();
		const volScalarField& alpha = thermo.alpha();

		//the heat flux of each surfaces of each volumes
		heatFlux = fvc::interpolate(alpha)*fvc::snGrad(he);

		//change the surface heat flux as the volume(vol) heat flux
		volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

		const surfaceScalarField::Boundary& heatFluxBf = heatFlux.boundaryField();

		//find the interface patch to calculate the heat flux
		wallHeatFluxBf[patchID] = heatFluxBf[patchID];
		wallHeatFlux = -wallHeatFlux;
		const scalarField& hfp = wallHeatFlux.boundaryField()[patchID];

		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			//Info<<hfp[faceI]<<endl;
		}

		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			//Find the location of preheat ergion
			if(y>pyrolysisUpperPosition){
				if ((hfp[faceI])<16000){
					Info<<"The location of preheat region is "<<y<<endl;
					preheatRegionUpperPosition = y;
					preheatLength = y-pyrolysisUpperPosition;
					if (mesh.time().value()<0.001){preheatLength = 0;}
					break;
				}
			}
		}

		currentStepTime = mesh.time().value();
		timeInterval = currentStepTime-previousStepTime;
		previousStepTime = currentStepTime;

		//construct condensed condensed region
		condensedRegion solidPanel(gMax(hfp), timeInterval, preheatLength);
		solidPanel.pyrolysisLengthIncrease_thermally_thin();
		if (mesh.time().value()>0.5)
		{
			pyrolysisUpperPosition = pyrolysisUpperPosition+solidPanel.pyrolysisLengthIncrease_;
		}
		solidPanel.pyrolysisBurnout();
		if (mesh.time().value()>3)
		{
			pyrolysisLowerPosition = pyrolysisLowerPosition+solidPanel.pysolysisBurnOutIncrease_;
		}

		//pyrolysisUpperPosition = 0.13;

		std::ofstream file_pyrolysisUpperPosition;

		file_pyrolysisUpperPosition.open("pyrolysisUpperPosition", std::ofstream::out | std::ofstream::app);

		file_pyrolysisUpperPosition<<runTime.timeName()<<token::TAB<<\
		preheatLength<<token::TAB\
		<<timeInterval<<token::TAB\
		<<pyrolysisUpperPosition<<token::TAB\
		<<solidPanel.pyrolysisLengthIncrease_<<token::TAB\
		<<gMax(hfp)<<token::TAB
		<<pyrolysisLowerPosition<<token::TAB<<std::endl;

		file_pyrolysisUpperPosition.close();

		std::ofstream file_pyrolysis;

		file_pyrolysis.open("pyrolysis", std::ofstream::out | std::ofstream::app);

		file_pyrolysis<<pyrolysisUpperPosition-0.13<<std::endl;

		file_pyrolysis.close();

		//update the fuel
		//Define the fuel boundary condition object
		fvPatchScalarField& fuel_patch = refCast<fvPatchScalarField>(composition.Y("C7H16").boundaryFieldRef()[patchID]);

		//Cast calculated fuel to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.002;
			scalar y_original = y-pyrolysisLowerPosition;

			if (y<pyrolysisLowerPosition){
				fuel_patch[faceI] = scalar(0);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Get the fuel
				nearWallRegion.calculateFuel();
				scalar meanInterFuel = nearWallRegion.interfaceFuel_;
				meanInterFuel = meanInterFuel + 0.1*meanInterFuel*Foam::sin(mesh.time().value());
				Random perturbation(rand());
				scalar interfaceFuel = meanInterFuel + 0.1*meanInterFuel*perturbation.GaussNormal();
				if(interfaceFuel<0){interfaceFuel = 0;}
				//Set fuel
				fuel_patch[faceI] = scalar( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition) ? meanInterFuel : 0 );
			}
		}
		Info<<"The fuel changes"<<endl;
		Info<<pyrolysisUpperPosition<<endl;


		//Define the velocity boundary condition object
		fvPatchVectorField& U_patch = refCast<fvPatchVectorField>(U.boundaryFieldRef()[patchID]);

		//Cast calculated velocity to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;
			//Info<<"Test"<<endl;

			if (y<pyrolysisLowerPosition){
				U_patch[faceI] = vector(0, 0, 0);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Info<<"value is "<<nearWallRegion.testValue<<endl;
				//Get the velocity
				nearWallRegion.calculateStreamwiseVelocity();
				scalar meanstreamwiseVelocity = nearWallRegion.streamwiseVelocity_;
				meanstreamwiseVelocity = meanstreamwiseVelocity + 0.1*meanstreamwiseVelocity*Foam::sin(mesh.time().value());
				Random perturbation_2(rand());
				scalar streamwiseVelocity = meanstreamwiseVelocity + 0.1*meanstreamwiseVelocity*perturbation_2.GaussNormal();
				//Info<<"The streamwiseVelocity is "<<streamwiseVelocity<<endl;
				nearWallRegion.calculateTransverseVelocity();
				scalar meantransverseVelocity = nearWallRegion.transverseVelocity_;
				meantransverseVelocity = meantransverseVelocity + 0.1*meantransverseVelocity*Foam::sin(mesh.time().value());
				Random perturbation_3(rand());
				scalar transverseVelocity = meantransverseVelocity + 0.1*meantransverseVelocity*perturbation_3.GaussNormal();
				Random perturbation_4(rand());
				scalar background = 0.1*perturbation_4.GaussNormal();
				//Info<<"The transverseVelocity is "<<transverseVelocity<<endl;
				//Set the velocity
				U_patch[faceI] = ( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition)\
				 ? vector(meanstreamwiseVelocity, meantransverseVelocity, 0) : vector(0, 0, 0));
			}

		}
		Info<<"The velocity changes"<<endl;

		//update the temperature
		//Define the temperature boundary condition object
		/*
		fvPatchScalarField& T_patch = refCast<fvPatchScalarField>(T.boundaryFieldRef()[patchID]);

		//Cast calculated temperature to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;

			if (y<pyrolysisLowerPosition){
				T_patch[faceI] = scalar(293);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Get the temperature
				nearWallRegion.calculateTemperature();
				scalar interfaceTemperature = nearWallRegion.interfaceTemperature_;
				//Set the temperature
				T_patch[faceI] = scalar( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition) ? interfaceTemperature : 293 );
			}

		}
		Info<<"The temperature changes"<<endl;
		*/

		//update the oxygen
		//Define the oxygen boundary condition object
		fvPatchScalarField& oxygen_patch = refCast<fvPatchScalarField>(composition.Y("O2").boundaryFieldRef()[patchID]);

		//Cast calculated oxygen to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;

			if (y<pyrolysisLowerPosition){
				oxygen_patch[faceI] = scalar(0.23);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Get the oxygen
				nearWallRegion.calculateOxygen();
				scalar interfaceOxygen = nearWallRegion.interfaceOxygen_;
				//Set oxygen
				oxygen_patch[faceI] = scalar( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition) ? interfaceOxygen : 0.23 );
			}
		}
		Info<<"The oxygen changes"<<endl;

		if ((mesh.time().value()>200) && (mesh.time().value())<445){condensedRegion::C  = 30;}
		if (mesh.time().value()>445){condensedRegion::C  = VGREAT;}

		Info<<"The index is "<<	assimilationTimeIndex <<endl;
		Info<<"The time is "<<	data[assimilationTimeIndex][0] <<endl;
		Info<<"The observation is "<<data[assimilationTimeIndex][1] <<endl;
