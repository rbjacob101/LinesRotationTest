using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class MasterScript : MonoBehaviour {

    public enum Point { E1P1, E1P2, E2P1, E2P2 };

    public Vector2 E1p1, E1p2, E2p1, E2p2, e1p1_last, e1p2_last;
    public Point point1;
    public Point point2;

    //for testing line intersections
    public Point intPoint1;
    public Point intPoint2;
    public Point radialPoint;
    public float angle;

    public float arcSize;
    public bool endpointLabels;
    public int direction;

    public Vector2 P1
    { get { return P2V(point1); } }
    public Vector2 P2
    { get { return P2V(point2); } }

    public Vector2 P2V(Point p)
    {
        switch (p)
        {
            case Point.E1P1:
                return E1p1;
            case Point.E1P2:
                return E1p2;
            case Point.E2P1:
                return E2p1;
            case Point.E2P2:
                return E2p2;
            default:
                return Vector2.zero;
        }
    }

    //calculates the angle between two edges going in a direction
    //dir = -1 -> ccw, dir = 1 -> cw
    public float EdgeEdgeAngle(Vector2 p11, Vector2 p12, Vector2 p21, Vector2 p22, int dir)
    {
        //if e1 and e2 are identical, return max angle 360
        if ((p11 == p21 && p12 == p22) || (p11 == p22 && p12 == p21))
            return 360;

        //raw square magnitudes, e.g. s11 = |p11|^2
        //calculate these only once
        float s11 = p11.sqrMagnitude, s12 = p12.sqrMagnitude,
            s21 = p21.sqrMagnitude, s22 = p22.sqrMagnitude;

        //vector2s in order of magnitude
        //|s11| <= |s12|, |s21| <= |s22|
        Vector2 e11, e12, e21, e22;
        if (s11 < s12) { e11 = p11; e12 = p12; }
        else { e11 = p12; e12 = p11; }
        if (s21 < s22) { e21 = p21; e22 = p22; }
        else { e21 = p22; e22 = p21; }

        //reorder the square magnitudes
        //se11 <= se12, se21 <= se22
        float se11 = Mathf.Min(s11, s12), se12 = Mathf.Max(s11, s12),
        se21 = Mathf.Min(s21, s22), se22 = Mathf.Max(s21, s22);

        //test if both points on e2 are outside of the annulus formed by e1's rotation
        //Approximately used for floating point errors
        if (((se21 > se12 || Mathf.Approximately(se21, se12)) && (se22 > se12 || Mathf.Approximately(se22, se12))) ||
            ((se21 < se11 || Mathf.Approximately(se21, se11)) && (se22 < se11 || Mathf.Approximately(se22, se11))))
            return 360;

        /* Test for points that share an endpoint endpoint
         * orb = 0 => no endpoint shared
         * orb = 1 => p1 endpoint orbitals (e12 and e22 shared)
         * orb = 2 => p2 endpoint orbitals (e11 and e21 shared) */
        sbyte orb = 0;
        if (e12 == e22) orb = 1;
        if (e11 == e21) orb = 2;

        //if there are orbitals, resolve their angles
        if (orb > 0)
        {
            //define the orbitals, o1 ∈ E1 and o2 ∈ E2
            Vector2 o1, o2;
            if (orb == 1) { o1 = e11; o2 = e21; }
            else { o1 = e12; o2 = e22; }

            //sign of determinant of o1, o2 determines clockwiseness of angle
            //1 = ccw, -1 = cw
            sbyte s = (sbyte)Mathf.Sign(o1.x * o2.y - o1.y * o2.x);

            //if sign of det multiplied by direction less than zero, rotating causes overlap, so angle = 0
            //else, more calculations needed
            if (dir * s < 0)
            {
                //if the shared endpoint is at (0,0),
                //angle is the angle between the orbitals
                if ((orb == 1 && e12 == Vector2.zero) ||
                    (orb == 2 && e11 == Vector2.zero))
                    return Vector2.Angle(o1, o2);
                return 0;
            }
            else
            {
                //calculate the square magnitudes of the orbitals for the purposes of comparison
                float o1s = o1.sqrMagnitude, o2s = o2.sqrMagnitude;

                //if |o1| equals |o2|, no need to calculate an intersection
                if (o1s == o2s)
                    return 360 - Vector2.Angle(o1, o2);

                //if |o1| and |o2| unequal, need to calculate an intersection
                if (o1s > o2s)
                    return 360 - Vector2.Angle(o2, EdgeCircleIntersection(o2.magnitude, e11, e12));
                return 360 - Vector2.Angle(o1, EdgeCircleIntersection(o1.magnitude, e21, e22));
            }
        }

        //independent lines, so need to find two angles and compare
        //near (n1, n2) and far (f1, f2) lines
        //n1, f1 ∈ E1 and n2, f2 ∈ E2
        Vector2 n1, n2, f1, f2;

        //near points
        if (se11 == se21) { n1 = e11; n2 = e21; }
        else if (se11 > se21) { n1 = e11; n2 = EdgeCircleIntersection(e11.magnitude, e21, e22); }
        else { n1 = EdgeCircleIntersection(e21.magnitude, e11, e12); n2 = e21; }

        //far points
        if (se12 == se22) { f1 = e12; f2 = e22; }
        else if (se12 > se22) { f1 = EdgeCircleIntersection(e22.magnitude, e11, e12); f2 = e22; }
        else { f1 = e12; f2 = EdgeCircleIntersection(e12.magnitude, e21, e22); }

        //sign of determinant of n1, n2 determines clockwiseness of angle
        //1 = ccw, -1 = cw
        //but if n1, n2 determinant = 0, check sign of determinant of f1, f2 as backup
        //as det = 0 does not give enough information about direction of angle
        float nDet = n1.x * n2.y - n1.y * n2.x;
        sbyte si = (sbyte)Mathf.Sign(nDet != 0 ? nDet : f1.x * f2.y - f1.y * f2.x);

        //determines whether the near angle (∠n1n2) is smaller
        bool nSmaller = false;

        //dot products of angles' rays (d1 corresponds to near angle, d2 to far)
        float d1 = Vector2.Dot(n1, n2), d2 = Vector2.Dot(f1, f2);
        if (d1 >= 0 && d2 <= 0)
            nSmaller = true;
        else
        {
            //floats used to compare cos^2 values of angles
            //because (d1*|f1|*|f2|)^2 = (|n1|*|n2|*|f1|*|f2|*cos ϴ)^2,
            //cos^2 ϴ values can be compared for dot products of same sign
            float rn = d1 * d1 * f1.sqrMagnitude * f2.sqrMagnitude,
                rf = d2 * d2 * n1.sqrMagnitude * n2.sqrMagnitude;

            //only in these two cases is the n angle (∠n1n2) smaller
            if ((d1 > 0 && d2 > 0 && rn > rf) ||
                (d1 < 0 && d2 < 0 && rn < rf))
                nSmaller = true;
        }

        //if sign of det multiplied by direction less than zero, rotating causes overlap, so use the smallest angle
        if (dir * si < 0)
            return nSmaller ? Vector2.Angle(n1, n2) : Vector2.Angle(f1, f2);
        //rotating away from the edge, so do 360 minus the bigger angle
        else return nSmaller ? 360 - Vector2.Angle(f1, f2) : 360 - Vector2.Angle(n1, n2);
    }

    //calculates the intersection between a circle centered around (0,0) with radius = r
    //and a line of form y = r or x = r, defined by two points p1 and p2
    public Vector2 EdgeCircleIntersection (float r, Vector2 p1, Vector2 p2)
    {
        if (p1.x == p2.x) {
            //vux
            float y = Mathf.Sqrt(r * r - p1.x * p1.x);
            if ((y <= p1.y || y <= p2.y) && (y >= p1.y || y >= p2.y)) return new Vector2(p1.x, y);
            else if ((-y <= p1.y || -y <= p2.y) && (-y >= p1.y || -y >= p2.y)) return new Vector2(p1.x, -y);
        } else {
            //hoy
            float x = Mathf.Sqrt(r * r - p1.y * p1.y);
            if ((x <= p1.x || x <= p2.x) && (x >= p1.x || x >= p2.x)) return new Vector2(x, p1.y);
            else if ((-x <= p1.x || -x <= p2.x) && (-x >= p1.x || -x >= p2.x)) return new Vector2(-x, p1.y);
        }
        return Vector2.zero;
    }
}

[CustomEditor(typeof(MasterScript))]
public class MasterScriptEditor : Editor
{
    public override void OnInspectorGUI ()
    {
        MasterScript s = target as MasterScript;
        serializedObject.Update();

        if (GUILayout.Button("Round to Integer"))
        {
            s.E1p1 = new Vector2(Mathf.RoundToInt(s.E1p1.x), Mathf.RoundToInt(s.E1p1.y));
            s.E1p2 = new Vector2(Mathf.RoundToInt(s.E1p2.x), Mathf.RoundToInt(s.E1p2.y));
            s.E2p1 = new Vector2(Mathf.RoundToInt(s.E2p1.x), Mathf.RoundToInt(s.E2p1.y));
            s.E2p2 = new Vector2(Mathf.RoundToInt(s.E2p2.x), Mathf.RoundToInt(s.E2p2.y));
            SceneView.RepaintAll();
        }

        EditorGUILayout.PropertyField(serializedObject.FindProperty("arcSize"), new GUIContent("Angle Arc Size"));
        EditorGUILayout.PropertyField(serializedObject.FindProperty("endpointLabels"), new GUIContent("Endpoint Labels?"));
        GUILayout.Space(10);
        EditorGUILayout.PropertyField(serializedObject.FindProperty("E1p1"), new GUIContent("E1P1"));
        EditorGUILayout.PropertyField(serializedObject.FindProperty("E1p2"), new GUIContent("E1P2"));
        EditorGUILayout.PropertyField(serializedObject.FindProperty("E2p1"), new GUIContent("E2P1"));
        EditorGUILayout.PropertyField(serializedObject.FindProperty("E2p2"), new GUIContent("E2P2"));
        GUILayout.Space(10);
        EditorGUILayout.BeginHorizontal();
        EditorGUIUtility.labelWidth = 100f;
        EditorGUILayout.PropertyField(serializedObject.FindProperty("point1"), new GUIContent("Rot. Start"));
        EditorGUILayout.PropertyField(serializedObject.FindProperty("point2"), new GUIContent("Rot. End"));
        EditorGUILayout.EndHorizontal();

        EditorGUILayout.BeginHorizontal();
        if (GUILayout.Button("Reset"))
        {
            s.E1p1 = s.e1p1_last;
            s.E1p2 = s.e1p2_last;
            SceneView.RepaintAll();
        }
        if (GUILayout.Button("Rotate"))
        {
            s.e1p1_last = s.E1p1;
            s.e1p2_last = s.E1p2;

            var a = Vector2.SignedAngle(s.P1, s.P2) * Mathf.Deg2Rad;
            s.E1p1 = new Vector2((s.E1p1.x * Mathf.Cos(a)) - (s.E1p1.y * Mathf.Sin(a)), (s.E1p1.y * Mathf.Cos(a)) + (s.E1p1.x * Mathf.Sin(a))); //rotate the points
            s.E1p2 = new Vector2((s.E1p2.x * Mathf.Cos(a)) - (s.E1p2.y * Mathf.Sin(a)), (s.E1p2.y * Mathf.Cos(a)) + (s.E1p2.x * Mathf.Sin(a)));
            SceneView.RepaintAll();
        }
        EditorGUILayout.EndHorizontal();

        EditorGUILayout.BeginHorizontal();
        if (GUILayout.Button("-2°"))
        {
            var a = (-2) * Mathf.Deg2Rad;
            s.E1p1 = new Vector2((s.E1p1.x * Mathf.Cos(a)) - (s.E1p1.y * Mathf.Sin(a)), (s.E1p1.y * Mathf.Cos(a)) + (s.E1p1.x * Mathf.Sin(a))); //rotate the points
            s.E1p2 = new Vector2((s.E1p2.x * Mathf.Cos(a)) - (s.E1p2.y * Mathf.Sin(a)), (s.E1p2.y * Mathf.Cos(a)) + (s.E1p2.x * Mathf.Sin(a)));
            SceneView.RepaintAll();
        }
        if (GUILayout.Button("+2°"))
        {
            var a = (2) * Mathf.Deg2Rad;
            s.E1p1 = new Vector2((s.E1p1.x * Mathf.Cos(a)) - (s.E1p1.y * Mathf.Sin(a)), (s.E1p1.y * Mathf.Cos(a)) + (s.E1p1.x * Mathf.Sin(a))); //rotate the points
            s.E1p2 = new Vector2((s.E1p2.x * Mathf.Cos(a)) - (s.E1p2.y * Mathf.Sin(a)), (s.E1p2.y * Mathf.Cos(a)) + (s.E1p2.x * Mathf.Sin(a)));
            SceneView.RepaintAll();
        }

        EditorGUILayout.EndHorizontal();

        EditorGUILayout.BeginHorizontal();
        if (GUILayout.Button("-5°"))
        {
            var a = (-5) * Mathf.Deg2Rad;
            s.E1p1 = new Vector2((s.E1p1.x * Mathf.Cos(a)) - (s.E1p1.y * Mathf.Sin(a)), (s.E1p1.y * Mathf.Cos(a)) + (s.E1p1.x * Mathf.Sin(a))); //rotate the points
            s.E1p2 = new Vector2((s.E1p2.x * Mathf.Cos(a)) - (s.E1p2.y * Mathf.Sin(a)), (s.E1p2.y * Mathf.Cos(a)) + (s.E1p2.x * Mathf.Sin(a)));
            SceneView.RepaintAll();
        }
        if (GUILayout.Button("+5°"))
        {
            var a = (5) * Mathf.Deg2Rad;
            s.E1p1 = new Vector2((s.E1p1.x * Mathf.Cos(a)) - (s.E1p1.y * Mathf.Sin(a)), (s.E1p1.y * Mathf.Cos(a)) + (s.E1p1.x * Mathf.Sin(a))); //rotate the points
            s.E1p2 = new Vector2((s.E1p2.x * Mathf.Cos(a)) - (s.E1p2.y * Mathf.Sin(a)), (s.E1p2.y * Mathf.Cos(a)) + (s.E1p2.x * Mathf.Sin(a)));
            SceneView.RepaintAll();
        }

        EditorGUILayout.EndHorizontal();

        EditorGUILayout.PropertyField(serializedObject.FindProperty("angle"), new GUIContent("Angle"));
        if (GUILayout.Button("Rotate")) {
            s.e1p1_last = s.E1p1;
            s.e1p2_last = s.E1p2;

            var a = (s.angle) * Mathf.Deg2Rad;
            s.E1p1 = new Vector2((s.E1p1.x * Mathf.Cos(a)) - (s.E1p1.y * Mathf.Sin(a)), (s.E1p1.y * Mathf.Cos(a)) + (s.E1p1.x * Mathf.Sin(a))); //rotate the points
            s.E1p2 = new Vector2((s.E1p2.x * Mathf.Cos(a)) - (s.E1p2.y * Mathf.Sin(a)), (s.E1p2.y * Mathf.Cos(a)) + (s.E1p2.x * Mathf.Sin(a)));
            SceneView.RepaintAll();
        }

        GUILayout.Space(10);

        EditorGUILayout.BeginHorizontal();
        if (GUILayout.Button("←"))
            s.direction = -1;
        if (GUILayout.Button("→"))
            s.direction = 1;
        EditorGUILayout.EndHorizontal();

        if (GUILayout.Button("Intersect Function"))
        {
            Debug.Log(s.EdgeEdgeAngle(s.E1p1, s.E1p2, s.E2p1, s.E2p2, s.direction));
        }

        GUILayout.Space(10);

        EditorGUILayout.BeginHorizontal();

        EditorGUILayout.PropertyField(serializedObject.FindProperty("intPoint1"), new GUIContent("1st Endpoint"));
        EditorGUILayout.PropertyField(serializedObject.FindProperty("intPoint2"), new GUIContent("2nd Endpoint"));

        EditorGUILayout.EndHorizontal();

        EditorGUILayout.PropertyField(serializedObject.FindProperty("radialPoint"), new GUIContent("Radial Point"));

        if (GUILayout.Button("Line Intersect Function"))
        {
            Vector2 point = s.EdgeCircleIntersection(s.P2V(s.radialPoint).magnitude, s.P2V(s.intPoint1), s.P2V(s.intPoint2));
            Debug.Log(string.Format("Intersection Point :: {0}", point));
        }

            serializedObject.ApplyModifiedProperties();
    }

    private void OnSceneGUI()
    {
        const float lineThickness = 3f;
        const float pointThickness = 0.02f;
        const float angleThreshold = 3f; //angle between P1 and P2 for arc graphic to show
        MasterScript s = target as MasterScript;

        Handles.color = new Color(1, 0, 0, 0.25f);

        //triangular sectors
        Handles.DrawAAConvexPolygon(Vector3.zero, s.E1p1, s.E1p2);
        Handles.DrawAAConvexPolygon(Vector3.zero, s.E2p1, s.E2p2);

        //labels
        Handles.Label((s.E1p1 + s.E1p2) / 2f, "L1");
        Handles.Label((s.E2p1 + s.E2p2) / 2f, "L2");
        Handles.Label(Vector3.zero, "C");
        if (s.endpointLabels)
        {
            Handles.Label(s.E1p1, "e1p1");
            Handles.Label(s.E1p2, "e1p2");
            Handles.Label(s.E2p1, "e2p1");
            Handles.Label(s.E2p2, "e2p2");
        }

        Handles.color = Color.white;

        //lines
        Handles.DrawAAPolyLine(lineThickness, s.E1p1, s.E1p2);
        Handles.DrawAAPolyLine(lineThickness, s.E2p1, s.E2p2);

        //points
        Handles.DrawSolidDisc(s.E1p1, Vector3.forward, pointThickness);
        Handles.DrawSolidDisc(s.E1p2, Vector3.forward, pointThickness);
        Handles.DrawSolidDisc(s.E2p1, Vector3.forward, pointThickness);
        Handles.DrawSolidDisc(s.E2p2, Vector3.forward, pointThickness);
        Handles.DrawSolidDisc(Vector3.zero, Vector3.forward, pointThickness);

        //annulus
        Handles.color = Color.red;
        Handles.DrawWireDisc(Vector3.zero, Vector3.forward, s.E1p1.magnitude);
        Handles.DrawWireDisc(Vector3.zero, Vector3.forward, s.E1p2.magnitude);

        //arc
        if (Vector2.Angle(s.P1, s.P2) >= angleThreshold)
        {
            Handles.color = Color.white;
            Handles.DrawAAPolyLine(lineThickness, s.P1.normalized * s.arcSize, Vector2.zero);
            Handles.DrawAAPolyLine(lineThickness, s.P2.normalized * s.arcSize, Vector2.zero);
            Handles.DrawWireArc(Vector3.zero, Vector3.forward, s.P1, Vector2.SignedAngle(s.P1, s.P2), s.arcSize);
            Handles.DrawWireArc(Vector3.zero, Vector3.forward, s.P1, Vector2.SignedAngle(s.P1, s.P2), s.arcSize - 0.005f);
            Handles.DrawWireArc(Vector3.zero, Vector3.forward, s.P1, Vector2.SignedAngle(s.P1, s.P2), s.arcSize - 0.01f);
            float om = ((s.P1.normalized * s.arcSize) - (s.P2.normalized * s.arcSize)).magnitude;
            float p1m = (s.P1.normalized * s.arcSize).magnitude;
            float p2m = (s.P2.normalized * s.arcSize).magnitude;
            Handles.Label(new Vector2(((p1m * (s.P2.normalized * s.arcSize).x) + (p2m * (s.P1.normalized * s.arcSize).x)) / (om + p1m + p2m),
                ((p1m * (s.P2.normalized * s.arcSize).y) + (p2m * (s.P1.normalized * s.arcSize).y)) / (om + p1m + p2m)),
                new GUIContent("θ"));
        }

    }
}
