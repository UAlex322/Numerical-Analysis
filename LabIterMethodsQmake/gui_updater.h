#pragma once
#include <QObject>

class GUIUpdater : public QObject {
    Q_OBJECT

public:
    explicit GUIUpdater(QObject *parent = 0) : QObject(parent) {}
signals:
    void requestUpdateCycles(const QString &);
    void finished();
};
