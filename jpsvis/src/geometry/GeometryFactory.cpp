#include "GeometryFactory.h"

#include "FacilityGeometry.h"
#include <vtkRenderer.h>
#include <vtkAssembly.h>

GeometryFactory::GeometryFactory()
{

}

void GeometryFactory::Init(vtkRenderer* renderer)
{
    for (auto&& rooms: _geometryFactory)
    {
        for(auto&& subroom:rooms.second)
        {
            subroom.second->CreateActors();
            renderer->AddActor(subroom.second->getActor2D());
            renderer->AddActor(subroom.second->getActor3D());
        }
    }
}

void GeometryFactory::Set2D(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->set2D(status);
        }
    }
    //Set3D(!status);
}

void GeometryFactory::Set3D(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->set3D(status);
        }
    }
    //Set2D(!status);
}

void GeometryFactory::Clear()
{
    for (auto&& src: _geometryFactory)
    {
        cout<<"cleaning...";
    }
}

void GeometryFactory::ChangeWallsColor(double* color)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->changeWallsColor(color);
        }
    }
}

void GeometryFactory::ChangeExitsColor(double* color)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->changeExitsColor(color);
        }
    }
}

void GeometryFactory::ChangeNavLinesColor(double* color)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->changeNavLinesColor(color);
        }
    }
}

void GeometryFactory::ChangeFloorColor(double* color)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->changeFloorColor(color);
        }
    }
}

void GeometryFactory::ChangeObstaclesColor(double* color)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->changeObstaclesColor(color);
        }
    }
}

void GeometryFactory::ShowDoors(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->showDoors(status);
        }
    }
}

void GeometryFactory::ShowStairs(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->showStairs(status);
        }
    }
}

void GeometryFactory::ShowWalls(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->showWalls(status);
        }
    }
}

void GeometryFactory::ShowNavLines(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->showNavLines(status);
        }
    }
}

void GeometryFactory::ShowFloor(bool status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->showFloor(status);
        }
    }
}

void GeometryFactory::ShowGeometryLabels(int status)
{
    for (auto&& room: _geometryFactory)
    {
        for(auto&& subroom:room.second)
        {
            subroom.second->showGeometryLabels(status);
        }
    }
}

bool GeometryFactory::RefreshView()
{
    if(_model.objectName()!="initialized")
    {
        _model.setObjectName("initialized");
        _model.setHorizontalHeaderItem( 0, new QStandardItem( "Entity" ) );
        //_model.setHorizontalHeaderItem( 1, new QStandardItem( "Description" ) );

        for (auto&& room: _geometryFactory)
        {
            //room caption
            //QStandardItem *roomcaption = new QStandardItem( QString("R %0").arg(room.first));
            //roomcaption->setEditable( false );
            //_model.setItem(room.first, 1, roomcaption);

            QStandardItem *item = new QStandardItem( QString("Room:%0").arg(room.first));
            item->setCheckable(true);
            item->setCheckState(Qt::Checked);

            for(auto&& subroom:room.second)
            {
                QStandardItem *child = new QStandardItem(
                            QString("%0: (%1)")
                            .arg(QString::fromStdString(subroom.second->GetDescription()))
                            .arg(subroom.first)
                            );

                child->setEditable(false);
                child->setCheckable(true);
                child->setCheckState(Qt::Checked);
                item->appendRow( child );
                _model.setItem(room.first, 0, item);
                QString data = QString("%0:%1").arg(room.first).arg(subroom.first);
                child->setData(data);

                //Subroom caption
                //QStandardItem *childcaption = new QStandardItem( QString("S %0").arg(subroom.first));
                //childcaption->setEditable( false );
                //_model.setItem(room.first, 1, childcaption);
            }
        }
        return true;
    }
    return false;
}

const std::map<int , std::map<int, std::shared_ptr<FacilityGeometry> > > & GeometryFactory::GetGeometry() const
{
    return _geometryFactory;
}

void GeometryFactory::AddElement(int room, int subroom, std::shared_ptr<FacilityGeometry> geo)
{
    //_geometryFactory.insert({room,{subroom,geo}});
    _geometryFactory[room][subroom]=geo;
}

std::shared_ptr<FacilityGeometry> GeometryFactory::GetElement(int room, int subroom)
{
    if(_geometryFactory.count(room))
    {
        if(_geometryFactory[room].count(subroom))
        {
            return _geometryFactory[room][subroom];
        }
    }
    return nullptr;
}

void GeometryFactory::UpdateVisibility(int room,int subroom,bool status)
{
    if(_geometryFactory.count(room))
    {
        if(_geometryFactory[room].count(subroom))
        {
            _geometryFactory[room][subroom]->setVisibility(status);
        }
    }
}

QStandardItemModel& GeometryFactory::GetModel()
{
    return _model;
}
