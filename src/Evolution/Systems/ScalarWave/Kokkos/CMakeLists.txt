#Distributed under the MIT License.
#See LICENSE.txt for details.

set(LIBRARY ScalarWaveKokkos)

    add_spectre_library(${LIBRARY})

        spectre_target_sources(
            ${LIBRARY} PRIVATE ApplyBoundaryCorrectionsToTimeDerivativeKokkos
                .cpp Batched /
            ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched
                .cpp Batched /
            ComputeInternalBoundaryTermsBatched.cpp Batched /
            ComputeTimeDerivativeBatched.cpp Batched /
            LiftInternalBoundaryTermsBatched.cpp Batched /
            ObserveTimeStepBatchedEvent.cpp Batched /
            PackageLocalFacesBatched.cpp ComputeTimeDerivativeKokkos
                .cpp RecordTimeStepperDataKokkos3D.cpp UpdateUKokkos3D.cpp)

            spectre_target_headers(
                ${LIBRARY} INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR} /
                src HEADERS
                    ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp Batched /
                Actions.hpp Batched /
                ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched
                    .hpp Batched /
                CleanHistoryBatched.hpp Batched /
                ComputeInternalBoundaryTermsBatched.hpp Batched /
                ComputeTimeDerivativeBatched.hpp Batched /
                DriverState.hpp Batched / InitializeBatchedData.hpp Batched /
                InitializeBoundaryBatchMetadata.hpp Batched /
                LiftInternalBoundaryTermsBatched.hpp Batched /
                ObserveTimeStepBatchedEvent.hpp Batched /
                PackageLocalFacesBatched.hpp Batched /
                RecordTimeStepperDataBatched.hpp Batched /
                UpdateUBatched.hpp CleanHistoryKokkos
                    .hpp ComputeTimeDerivativeKokkos
                    .hpp InitializeKokkosBoundaryCommunication
                    .hpp InitializeKokkosTags
                    .hpp InitializeKokkosTimeStepperState
                    .hpp KokkosBoundaryCommunication.hpp MortarDataPrimitives
                    .hpp KokkosTimeStepperTags.hpp RecordTimeStepperDataKokkos
                    .hpp RecordTimeStepperDataKokkos3D
                    .hpp SendBoundaryCorrectionsKokkos.hpp SyncKokkosToHost
                    .hpp UpdateUKokkos.hpp UpdateUKokkos3D.hpp)

                target_link_libraries(
                    ${LIBRARY} PUBLIC DataStructures Domain ErrorHandling
                        EventsAndTriggers Options Parallel ScalarWave
                            Serialization Time TimeKokkos Utilities PRIVATE
                                LinearOperators)
