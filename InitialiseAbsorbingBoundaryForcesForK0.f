        subroutine InitialiseAbsorbingBoundaryForcesForK0(NodalVerticalStressesSoil, NodalPressuresWater)
        !************************************************************************************
        !
        !   Function:  Calculate the dashpot terms
        !
        !************************************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%Nodtot), intent(in) :: NodalVerticalStressesSoil, NodalPressuresWater
          ! Local variables
          integer(INTEGER_TYPE) :: INod, NiX, IEntity, IDim
          real(REAL_TYPE) :: A, F, FX, FY, FZ, FWP, FWPX, FWPY, FWPZ, K0Value
          real(REAL_TYPE), dimension(NVECTOR):: IConditionXYZ
 
          K0Value = MatParams(CalParams%AbsorbingBoundaries%VBMaterialSet)%K0Value
      
          do INod = 1, Counters%NodTot
          
            A = VisNodSurAraSld(INod) ! Surface area
            if (A>0.0) then ! Get stresses at this node
              F = NodalVerticalStressesSoil(INod) * A
              FWP = NodalPressuresWater(INod) * A
              NiX = ReducedDof(INod) + 1
              
              do IDim = 1, NVECTOR !over coordinates
                IConditionXYZ(IDim) = 	NodeCndVisSld2(INod, IDim)
              end do

              if (IConditionXYZ(1)==1) then ! Normal direction ... give force 
                FX = F * K0Value
                FWPX = FWP
              else
                FX = 0.0
                FWPX = 0.0
              end if
        
              if (IConditionXYZ(2)==1) then ! Normal direction ... give force 
                FY = F
                FWPY = FWP
              else
                FY = 0.0
                FWPY = 0.0
              end if
               
              if (NVECTOR==3) then			   
                if (IConditionXYZ(3)==1) then ! Normal direction ... give force
                  FZ = F * K0Value
                  FWPZ = FWP
                else
                  FZ = 0.0
                  FWPZ = 0.0
                end if
              end if 
			  
              do IEntity = 1, Counters%NEntity
                VisDampForceSld(NiX, IEntity) = -FX
                VisDampForceSld(NiX + 1, IEntity) = FY
                VisDampForceWat(NiX, IEntity) = FWPX
                VisDampForceWat(NiX + 1, IEntity) = FWPY

                if (NVECTOR==3) then	
                  VisDampForceSld(NiX + 2, IEntity) = FZ		
                  VisDampForceWat(NiX + 2, IEntity) = FWPZ
                end if 			
	
              end do
                
            end if
          end do

          if (IS3DCYLINDRIC) then
           do IEntity = 1, Counters%nEntity
             call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, VisDampForceSld(:, IEntity), VisDampForceSld(:, IEntity))
           end do
         end if

        end subroutine InitialiseAbsorbingBoundaryForcesForK0