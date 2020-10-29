
#include "Temperature.H"
#include "primitivePatchInterpolation.H"
#include "volFields.H"

#include "apiCoupledTemperatureFvPatchScalarField.H"

using namespace Foam;

preciceAdapter::CHT::Temperature::Temperature
(
    const Foam::fvMesh& mesh,
    const std::string nameT
) :
    CouplingDataUser(DT_Scalar),
    T_(mesh.lookupObjectRef<volScalarField>(nameT)),
    mesh_(mesh)
{}

void preciceAdapter::CHT::Temperature::write(std::vector<double> &buffer, bool meshConnectivity, const unsigned int dim)
{
    std::size_t bufferIndex = 0;

    // For every boundary patch of the interface
    for (std::size_t j = 0; j < patchIDs_.size(); j++ )
    {
        const auto   patchID          (patchIDs_.at(j));
        const auto & boundaryPatch    (T_.boundaryField()[patchID]);

        //If we use the mesh connectivity, we interpolate from the centres to the nodes
        if(meshConnectivity)
        {
            //Create an Interpolation object at the boundary Field
            primitivePatchInterpolation patchInterpolator(mesh_.boundaryMesh()[patchID]);

            //Interpolate from centers to nodes
            const scalarField pointValue (patchInterpolator.faceToPointInterpolate(boundaryPatch));

            std::cout << "CHT::Temperature::write (" << j << ") ";
            forAll(pointValue, i)
            {
                buffer[bufferIndex++] = pointValue[i];
                std::cout << boundaryPatch[i] << ", ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cout << "CHT::Temperature::write (" << j << ") ";
            forAll(boundaryPatch, i)
            {
                buffer[bufferIndex++] = boundaryPatch[i];
                std::cout << boundaryPatch[i] << ", ";
            }
            std::cout << std::endl;
        }
    }
}

void preciceAdapter::CHT::Temperature::read(const std::vector<double> &buffer, const unsigned int dim)
{
    std::size_t bufferIndex = 0;

    // For every boundary patch of the interface
    for (std::size_t j = 0; j < patchIDs_.size(); j++)
    {
        const auto  patchID         (patchIDs_.at(j));
        auto &      boundaryPatch   (refCast<apiCoupledTemperatureFvPatchScalarField> (T_.boundaryFieldRef()[patchID]));
        auto &      value           (boundaryPatch.refValue());

        std::cout << "CHT::Temperature::read (" << j << ") ";
        forAll(value, i)
        {
            value[i] = buffer[bufferIndex++];
            std::cout << value[i] << ", ";
        }
        std::cout << std::endl;
    }
}
