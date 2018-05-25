#ifndef __37e385afbd3b4b1dba8611fb71787822
#define __37e385afbd3b4b1dba8611fb71787822
#include <lf/base/base.h>
#include "geometry.h"


namespace lf::mesh {

class Entity {

public:

  
  /**
   * @brief The codimension of this entity w.r.t. the Mesh.dimMesh()
   * of the owning mesh manager.
   */
  virtual char Codim() const = 0;

  /**
   * @brief Return all sub entities of this entity that have the given 
   *        codimension (w.r.t. this entity!)
   * @param codim The codim w.r.t. this entity
   * @return 
   */
  virtual base::RandomAccessRange<Entity> SubEntities(char codim) const = 0;

  /**
   * @brief Describes the geometry of this entity.
   * @return A point to a Geometry object that will remain valid for as long
   *         as the Mesh remains valid.
   */
  virtual Geometry* Geometry() const = 0;

  /**
   * @brief Describes the reference element type of this entity.
   * @return An object of type base::RefEl.
   */
  virtual base::RefEl RefEl() const = 0;

  virtual ~Entity() {}
};


}

#endif // __37e385afbd3b4b1dba8611fb71787822